// Linked against actual upstream objects; only the original CLI main is renamed.
#include "rlink.h"
#include <iomanip>
#include <iostream>
#include <string>

extern bool longreads, mixedMode, isunitig, nomulti;

void junction(std::ostream& out, const CJunction& j) {
    out << j.start << ',' << j.end << ',' << int(j.strand) << ',' << int(j.guide_match)
        << ',' << int(j.consleft) << ',' << int(j.consright) << ',' << j.nreads
        << ',' << j.nreads_good << ',' << j.leftsupport << ',' << j.rightsupport
        << ',' << j.nm << ',' << j.mm << ';';
}
void snapshot(const BundleData& b) {
    std::cout << std::hexfloat << "B:" << b.start << ',' << b.end << ',' << b.numreads
              << ',' << b.num_fragments << ',' << b.frag_len << ',' << b.sum_cov << ',' << int(b.covflags) << ',' << int(b.status);
    // GList/GVec accessors in this upstream version are not all const-qualified.
    BundleData& data=const_cast<BundleData&>(b);
    for (int i=0; i<data.readlist.Count(); ++i) {
        CReadAln& r=*data.readlist[i];
        std::cout << " R:" << i << ',' << r.start << ',' << r.end << ',' << int(r.strand)
                  << ',' << r.nh << ',' << r.len << ',' << r.read_count << ',' << r.unitig << ',' << r.longread << ',' << r.in_guide;
        for (int e=0; e<r.segs.Count(); ++e) std::cout << " E:" << r.segs[e].start << ',' << r.segs[e].end;
        for (int p=0; p<r.pair_idx.Count(); ++p) std::cout << " P:" << r.pair_idx[p] << ',' << r.pair_count[p];
        for (int j=0; j<r.juncs.Count(); ++j) { std::cout << " J:"; junction(std::cout, *r.juncs[j]); }
    }
    for (int j=0; j<data.junction.Count(); ++j) { std::cout << " G:"; junction(std::cout, *data.junction[j]); }
    for (int s=0; s<3; ++s) for (int i=0; i<data.bpcov[s].Count(); ++i) std::cout << " C:" << s << ',' << i << ',' << data.bpcov[s][i];
    std::cout << '\n';
}
int main(int argc, char** argv) {
    if (argc!=3) return 2;
    std::string mode=argv[2];
    longreads=mode=="long" || mode=="unitig";
    mixedMode=mode=="mixed";
    isunitig=mode=="unitig";
    nomulti=mode=="nomulti";
    mergeMode=false;
    BundleData bundle;
    GHash<int> pairs;
    GSamReader reader(argv[1]);
    GSamRecord* record;
    int reference=-1, end=0, calls=0;
    while ((record=reader.next())!=NULL) {
        if (record->refId()!=reference || int(record->start)>end+50) {
            bundle.getReady(bundle.start, bundle.end);
            snapshot(bundle);
            bundle.Clear(); pairs.Clear();
            reference=record->refId(); end=0;
        }
        record->uval=record->tag_int("ZU");
        const char xs=record->spliceStrand();
        GReadAlnData data(record, xs=='+' ? 1 : xs=='-' ? -1 : 0,
                         record->tag_int("NH") ? record->tag_int("NH") : 1,
                         record->tag_int("HI"));
        // Deliberately simulate both outcomes of the real guide evaluator.
        for (int i=1; i<record->exons.Count(); ++i) {
            CJunction* j=new CJunction(record->exons[i-1].end+1, record->exons[i].start-1);
            j->guide_match=record->tag_int("ZG");
            data.juncs.Add(j);
        }
        int start=bundle.readlist.Count() ? bundle.start : int(record->start);
        if (end<int(record->end)) end=record->end;
        processRead(start, end, bundle, pairs, data);
        std::cout << "CALL:" << ++calls << ' ';
        snapshot(bundle);
        delete record;
    }
    bundle.getReady(bundle.start, bundle.end);
    snapshot(bundle);
    bundle.Clear(); snapshot(bundle);
}
