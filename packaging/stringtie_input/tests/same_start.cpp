#include "same_start_index.h"
#include <cassert>
#include <iostream>
#include <random>

struct Seg { unsigned start, end; };
struct Exons : std::vector<Seg> {
    int Count() const { return static_cast<int>(size()); }
};
struct Read {
    unsigned start;
    char strand;
    bool longread, unitig;
    Exons segs;
    unsigned stored_del;
};
struct Reads : std::vector<Read*> {
    int Count() const { return static_cast<int>(size()); }
};

bool predicate(const Read& a, const Read& b, bool deletion, bool unitig) {
    if (a.strand != b.strand || a.longread != b.longread ||
        (unitig && a.unitig != b.unitig) || a.segs.size() != b.segs.size()) return false;
    for (size_t i=0; i<a.segs.size(); ++i)
        if (a.segs[i].start != b.segs[i].start || a.segs[i].end != b.segs[i].end) return false;
    return !deletion || a.stored_del == b.stored_del;
}
int linear(const Reads& reads, const Read& q, bool deletion, bool unitig) {
    for (int i=reads.Count()-1; i>=0 && reads[i]->start==q.start; --i)
        if (predicate(*reads[i], q, deletion, unitig)) return i;
    return -1;
}
int indexed(StringTieSameStartIndex& index, Reads& reads, Read& q, bool deletion, bool unitig) {
    const std::vector<int>& candidates=index.find(reads, q.segs, q.strand, q.longread, q.start);
    for (auto i=candidates.rbegin(); i!=candidates.rend(); ++i)
        if (predicate(*reads[*i], q, deletion, unitig)) return *i;
    return -1;
}
Read read(unsigned start=100, unsigned end=200) {
    Read r; r.start=start; r.strand=1; r.longread=true; r.unitig=false;
    r.stored_del=0; r.segs.push_back({start, end}); return r;
}
int main() {
    StringTieSameStartIndex index;
    Read a=read(), b=read(), c=read(), q=read();
    b.stored_del=5; c.stored_del=7;
    Reads reads; reads.push_back(&a); reads.push_back(&b); reads.push_back(&c);
    assert(indexed(index, reads, q, true, false)==0); // must keep searching past mismatches
    assert(indexed(index, reads, q, false, false)==2); // newest, not first inserted
    c.unitig=true; // mutable field is deliberately NOT part of the key
    assert(indexed(index, reads, q, false, true)==1);
    a.strand=-1; // only test mutation after lifecycle reset; key fields are input-immutable
    index.clear();
    q.strand=-1;
    assert(indexed(index, reads, q, true, false)==0);
    q.strand=0; assert(indexed(index, reads, q, false, false)==-1);
    q.strand=1; q.longread=false; assert(indexed(index, reads, q, false, false)==-1);
    q=read(101); assert(indexed(index, reads, q, false, false)==-1);
    Read later=q; reads.push_back(&later); assert(indexed(index, reads, q, false, false)==3);
    q=read(100); // a previous start is no longer in the trailing run
    assert(indexed(index, reads, q, false, false)==-1);
    reads.clear(); index.clear(); // same count and coordinates, different bundle
    a=read(); reads.push_back(&a);
    assert(indexed(index, reads, q, false, false)==0);

    std::mt19937 rng(222);
    std::vector<Read> storage; storage.reserve(12000);
    reads.clear(); index.clear();
    unsigned start=100;
    for (int i=0; i<10000; ++i) {
        if (i%250==0) ++start;
        Read next=read(start, start+100+rng()%20);
        next.strand=static_cast<char>(static_cast<int>(rng()%3)-1);
        next.longread=rng()%2; next.unitig=rng()%2;
        next.stored_del=rng()%4; // stored junction signature, not original input signature
        if (rng()%2) next.segs.push_back({start+500, start+600+static_cast<unsigned>(rng()%10)});
        for (int d=0; d<2; ++d) for (int u=0; u<2; ++u)
            assert(indexed(index, reads, next, d, u)==linear(reads, next, d, u));
        if (linear(reads, next, true, true)<0) {
            storage.push_back(next); reads.push_back(&storage.back());
        }
        if (i%701==0) { reads.clear(); index.clear(); }
    }

    // An exact-exon bucket avoids scanning unrelated same-start alignments.
    reads.clear(); storage.clear(); index.clear();
    for (unsigned i=0; i<6000; ++i) { storage.push_back(read(1, 100+i)); reads.push_back(&storage.back()); }
    q=read(1, 100);
    const std::vector<int>& candidates=index.find(reads, q.segs, q.strand, q.longread, q.start);
    assert(candidates.size()==1 && candidates[0]==0);
    std::cout << "candidate regressions passed (40000 differential queries)\n";
}
