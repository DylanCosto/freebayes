//
// Multi-threaded calling (--threads), included by freebayes.cpp
//
// The targets are split into chunks which worker threads process with their
// own AlleleParser, and the main thread writes the records in order. The
// output is the same as with a single thread:
//
//  - Calling at a position depends on the parser's history (e.g. a haplotype
//    call makes it skip the positions it covers), so a chunk is processed
//    starting PAD bases before its nominal start and continuing PAD bases
//    after its end, recording the positions it evaluates and the records it
//    produces.
//  - Two neighbouring chunks are joined at the first position which both
//    evaluate in their overlap (at least MARGIN bases from the ends of the
//    processed ranges). From there on both follow the same path; this is
//    verified over the rest of the overlap (same positions evaluated, same
//    records). The records before the join come from the first chunk, those
//    after it from the second.
//  - If the two cannot be joined, the second chunk is processed again
//    starting from the start of the first chunk, which reproduces the
//    single-threaded path exactly.
//
// The chunks have about the same amount of data, estimated from the BAM
// indexes, as the time taken depends mostly on the number of reads, and
// regions such as collapsed repeats can have a hundred times the usual depth.
//

#include <htslib/sam.h>

namespace fbparallel {

static const long PAD = 1000;
static const long MARGIN = 300;
static const long WINDOW = 16384;           // resolution of the data estimate
static const long MAX_CHUNK = 1000000;
static const long MIN_CHUNK = 2000;
static const int CHUNKS_PER_THREAD = 16;

struct Chunk {
    string seq;
    long tleft, tright;         // enclosing target, 0-based inclusive
    long left, right;           // nominal range
    long runLeft, runRight;     // processed range
    double weight;              // estimated amount of data
    CallLog log;
    int status;                 // 0: pending, 1: done
};

// Process [runLeft, runRight] of a chunk with a new parser
static void runChunk(const Parameters& base, Chunk& chunk) {
    Parameters params = base;
    params.regions.clear();
    params.targets = "";
    params.outputFile = "";
    params.threads = 1;
    ostringstream region;
    // freebayes regions are 0-based, end-exclusive
    region << chunk.seq << ":" << chunk.runLeft << "-" << chunk.runRight + 1;
    params.regions.push_back(region.str());
    AlleleParser* parser = new AlleleParser(params);
    ostringstream unused;
    chunk.log = CallLog();
    runCaller(parser, unused, &chunk.log);
    delete parser;
}

// Opened BAM files with their indexes, for estimating the amount of data
struct IndexedBams {
    vector<samFile*> files;
    vector<sam_hdr_t*> headers;
    vector<hts_idx_t*> indexes;
    bool ok;

    IndexedBams(const vector<string>& bams) : ok(true) {
        for (vector<string>::const_iterator f = bams.begin(); f != bams.end(); ++f) {
            samFile* fp = sam_open(f->c_str(), "r");
            sam_hdr_t* hdr = fp ? sam_hdr_read(fp) : NULL;
            hts_idx_t* idx = hdr ? sam_index_load(fp, f->c_str()) : NULL;
            files.push_back(fp);
            headers.push_back(hdr);
            indexes.push_back(idx);
            if (!idx) ok = false;
        }
    }

    ~IndexedBams() {
        for (size_t i = 0; i < files.size(); ++i) {
            if (indexes[i]) hts_idx_destroy(indexes[i]);
            if (headers[i]) sam_hdr_destroy(headers[i]);
            if (files[i]) sam_close(files[i]);
        }
    }

    // Compressed bytes of the alignments overlapping [beg, end), summed over
    // the files
    double weight(const string& seq, long beg, long end) {
        double w = 0;
        for (size_t i = 0; i < files.size(); ++i) {
            int tid = sam_hdr_name2tid(headers[i], seq.c_str());
            if (tid < 0) continue;
            hts_itr_t* itr = sam_itr_queryi(indexes[i], tid, beg, end);
            if (!itr) continue;
            for (int j = 0; j < itr->n_off; ++j) {
                // virtual offsets: the compressed offset is in the upper 48 bits
                w += (double) ((itr->off[j].v >> 16) - (itr->off[j].u >> 16)) + 1;
            }
            hts_itr_destroy(itr);
        }
        return w;
    }
};

static void makeChunks(const BedTarget& t, IndexedBams& bams, double target, vector<Chunk>& chunks) {
    Chunk c;
    c.seq = t.seq;
    c.tleft = t.left;
    c.tright = t.right;
    c.status = 0;
    c.left = t.left;
    c.weight = 0;
    for (long w = t.left; w <= t.right; w += WINDOW) {
        long wright = min((long) t.right, w + WINDOW - 1);
        double weight = bams.weight(t.seq, w, wright + 1);
        if (weight > target) {
            // a dense window: split it evenly
            if (c.left < w) {
                c.right = w - 1;
                chunks.push_back(c);
            }
            long n = min((long) ceil(weight / target), (wright - w + 1) / MIN_CHUNK);
            n = max(n, 1L);
            long size = (wright - w + 1 + n - 1) / n;
            for (long l = w; l <= wright; l += size) {
                c.left = l;
                c.right = min(wright, l + size - 1);
                c.weight = weight * (c.right - c.left + 1) / (wright - w + 1);
                chunks.push_back(c);
            }
            c.left = wright + 1;
            c.weight = 0;
            continue;
        }
        if (c.left < w && (c.weight + weight > target || wright - c.left + 1 > MAX_CHUNK)) {
            c.right = w - 1;
            chunks.push_back(c);
            c.left = w;
            c.weight = 0;
        }
        c.weight += weight;
    }
    if (c.left <= t.right) {
        c.right = t.right;
        chunks.push_back(c);
    }
}

// Join chunks a and b (b follows a in the same target): returns the position
// from which the records of b are used, or -1 if they cannot be joined
static long join(const Chunk& a, const Chunk& b) {
    long lo = b.runLeft + MARGIN;
    long hi = a.runRight - MARGIN;
    if (lo > hi) return -1;
    const vector<long>& pa = a.log.positions;
    const vector<long>& pb = b.log.positions;
    vector<long>::const_iterator ia = lower_bound(pa.begin(), pa.end(), lo);
    vector<long>::const_iterator ib = lower_bound(pb.begin(), pb.end(), lo);
    // first position evaluated by both
    while (ia != pa.end() && ib != pb.end() && *ia <= hi && *ib <= hi && *ia != *ib) {
        if (*ia < *ib) ++ia; else ++ib;
    }
    if (ia == pa.end() || ib == pb.end() || *ia > hi || *ib > hi) return -1;
    long q = *ia;
    // verify that both evaluate the same positions and give the same records up to hi
    while (ia != pa.end() && *ia <= hi) {
        if (ib == pb.end() || *ib != *ia) return -1;
        ++ia; ++ib;
    }
    if (ib != pb.end() && *ib <= hi) return -1;
    vector<pair<long, string> >::const_iterator ra = a.log.records.begin(), rb = b.log.records.begin();
    while (ra != a.log.records.end() && ra->first < q) ++ra;
    while (rb != b.log.records.end() && rb->first < q) ++rb;
    while (ra != a.log.records.end() && ra->first <= hi) {
        if (rb == b.log.records.end() || *rb != *ra) return -1;
        ++ra; ++rb;
    }
    if (rb != b.log.records.end() && rb->first <= hi) return -1;
    return q;
}

} // namespace fbparallel

static bool canRunParallel(AlleleParser* parser) {
    Parameters& p = parser->parameters;
    if (p.threads < 2) return false;
    if (p.useStdin || p.gVCFout || p.limitCoverage > 0 || p.output != "vcf") return false;
    if (!p.variantPriorsFile.empty() || !p.haplotypeVariantFile.empty()) return false;
    if (!fbparallel::IndexedBams(p.bams).ok) {
        cerr << "freebayes: the input is not indexed, running single-threaded" << endl;
        return false;
    }
    return true;
}

static int runParallel(AlleleParser* mainParser, ostream& out) {
    using namespace fbparallel;
    const Parameters& params = mainParser->parameters;
    IndexedBams bams(params.bams);

    // the targets, or all the reference sequences with reads
    vector<BedTarget> targets = mainParser->targets;
    if (targets.empty()) {
        for (REFVEC::iterator r = mainParser->referenceSequences.begin(); r != mainParser->referenceSequences.end(); ++r) {
            if (r->REFLEN > 0 && bams.weight(r->REFNAME, 0, r->REFLEN) > 0)
                targets.push_back(BedTarget(r->REFNAME, 0, r->REFLEN - 1));
        }
    }

    double total = 0;
    for (vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t)
        total += bams.weight(t->seq, t->left, t->right + 1);
    double target = max(total / (params.threads * CHUNKS_PER_THREAD), 1.0);

    vector<Chunk> chunks;
    for (vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t)
        makeChunks(*t, bams, target, chunks);
    for (vector<Chunk>::iterator c = chunks.begin(); c != chunks.end(); ++c) {
        c->runLeft = max(c->tleft, c->left - PAD);
        c->runRight = min(c->tright, c->right + PAD);
    }

    // Processing order: the chunks with several times the usual amount of
    // data first, as they take longest, then the others in order (so that
    // their records can be written as they are done)
    vector<size_t> order;
    for (size_t k = 0; k < chunks.size(); ++k)
        if (chunks[k].weight > 2 * target) order.push_back(k);
    for (size_t k = 0; k < chunks.size(); ++k)
        if (chunks[k].weight <= 2 * target) order.push_back(k);

    mutex lock;
    condition_variable cond;
    size_t next = 0;
    vector<thread> workers;
    int nthreads = min((int) chunks.size(), params.threads);
    for (int i = 0; i < nthreads; ++i) {
        workers.push_back(thread([&]() {
            while (true) {
                size_t k;
                {
                    unique_lock<mutex> lk(lock);
                    if (next >= order.size()) break;
                    k = order[next++];
                }
                runChunk(params, chunks[k]);
                {
                    unique_lock<mutex> lk(lock);
                    chunks[k].status = 1;
                }
                cond.notify_all();
            }
        }));
    }

    const auto wait_for = [&](size_t k) {
        unique_lock<mutex> lk(lock);
        cond.wait(lk, [&]() { return chunks[k].status != 0; });
    };

    long from = LONG_MIN;   // records of the current chunk are used from this position
    for (size_t k = 0; k < chunks.size(); ++k) {
        wait_for(k);
        long to = LONG_MAX;
        bool sameTarget = k + 1 < chunks.size()
            && chunks[k + 1].seq == chunks[k].seq && chunks[k + 1].tleft == chunks[k].tleft;
        if (sameTarget) {
            wait_for(k + 1);
            to = join(chunks[k], chunks[k + 1]);
            if (to < 0) {
                // Process the next chunk again from the start of this one: it
                // then follows the single-threaded path
                Chunk& n = chunks[k + 1];
                n.runLeft = chunks[k].runLeft;
                runChunk(params, n);
                to = n.left;
            }
        }
        for (vector<pair<long, string> >::iterator r = chunks[k].log.records.begin(); r != chunks[k].log.records.end(); ++r) {
            if (r->first >= from && r->first < to) out << r->second;
        }
        chunks[k].log = CallLog();
        from = to;
    }

    for (size_t i = 0; i < workers.size(); ++i) workers[i].join();
    return 0;
}
