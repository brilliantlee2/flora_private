#ifndef FLORA_STRINGTIE_SAME_START_INDEX_H
#define FLORA_STRINGTIE_SAME_START_INDEX_H

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <vector>

// Only immutable input fields belong in this key. In particular, unitig and
// guide-adjusted deletion junctions must still be checked by processRead.
class StringTieSameStartIndex {
    typedef std::vector<unsigned> Key;
    struct Hash {
        std::size_t operator()(const Key& key) const {
            std::size_t h=0;
            for (unsigned v : key) h^=std::size_t(v)+0x9e3779b9U+(h<<6)+(h>>2);
            return h;
        }
    };
    typedef std::unordered_map<Key, std::vector<int>, Hash> Buckets;
    Buckets buckets_;
    const std::vector<int> empty_;
    int indexed_=0;
    unsigned start_=0;
    bool active_=false;

    template<class Exons>
    static Key key(Exons& exons, char strand, bool longread) {
        Key result;
        result.reserve(3+2*exons.Count());
        result.push_back(static_cast<unsigned char>(strand));
        result.push_back(longread);
        result.push_back(exons.Count());
        for (int i=0; i<exons.Count(); ++i) {
            result.push_back(exons[i].start);
            result.push_back(exons[i].end);
        }
        return result;
    }

public:
    static bool enabled() {
        // Immutable for the process lifetime; unknown modes fall back to linear.
        static const bool use_index=[]() {
            const char* mode=std::getenv("STRINGTIE_INPUT_LOOKUP");
            return mode==NULL || std::strcmp(mode, "indexed")==0;
        }();
        return use_index;
    }

    void clear() {
        Buckets().swap(buckets_); // release the entire run, including bucket capacity
        indexed_=0;
        active_=false;
    }

    template<class Reads, class Exons>
    const std::vector<int>& find(Reads& reads, Exons& exons, char strand,
                                 bool longread, unsigned start) {
        if (!active_ || start!=start_ || reads.Count()<indexed_) {
            clear();
            active_=true;
            start_=start;
            indexed_=reads.Count();
            // Exactly the contiguous suffix visited by upstream's reverse loop.
            while (indexed_>0 && reads[indexed_-1]->start==start) --indexed_;
        }
        if (reads.Count()>0 && reads[reads.Count()-1]->start!=start) return empty_;
        for (; indexed_<reads.Count(); ++indexed_) {
            auto* read=reads[indexed_];
            buckets_[key(read->segs, read->strand, read->longread)].push_back(indexed_);
        }
        auto found=buckets_.find(key(exons, strand, longread));
        return found==buckets_.end() ? empty_ : found->second;
    }
};

#endif
