#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <memory>//shared_ptr

// debug
#include <iostream>

template<typename T>
uint32_t ilog22(T x) {
	uint32_t r = 0;

	for (/**/; x; ++r) x >>= 1;

	return r;
}

struct djinn_ctx_buf_t {
    djinn_ctx_buf_t() : buffer(nullptr), pos(nullptr), capac(0) {}
    djinn_ctx_buf_t(size_t cap) : buffer(new uint8_t[cap]), pos(buffer), capac(cap) {}
    ~djinn_ctx_buf_t() { delete[] buffer; buffer = nullptr; }

	void reset() { pos = buffer; }
	size_t size() const { return(pos - buffer); }

    uint8_t* buffer;
    uint8_t* pos;
    size_t capac;
};

// Based on Subbotin Range Coder.
class CBasicRangeCoder {
public:
    CBasicRangeCoder() : low(0), range(0){}
	virtual ~CBasicRangeCoder() {}

	virtual void Start() =0;
	virtual void End() =0;

public:
	static constexpr uint32_t TopValue = 1 << 24;
	static constexpr uint32_t Mask32 = (uint32_t)-1;

    uint64_t low;
	uint32_t range;

    std::shared_ptr<djinn_ctx_buf_t> data;
};

class CRangeEncoder : public CBasicRangeCoder {
public:
	CRangeEncoder(){}

	void Start() override {
		low = 0;
		range = Mask32;
	}

	void EncodeFrequency(uint32_t symFreq, uint32_t cumFreq, uint32_t totalFreqSum) {
		assert(range > totalFreqSum);
		low   += cumFreq * (range /= totalFreqSum);
    	range *= symFreq;
		// range /= totalFreqSum;
		// low += range * cumFreq;
		// range *= symFreq;

		while (range < TopValue) {
			assert(range != 0);
			if ( (uint8_t)((low ^ (low + range)) >> 56) )
            	range = (((uint32_t)(low) | (TopValue - 1)) - (uint32_t)(low));

            *data->pos++ = (uint8_t)(low >> 56);
			low <<= 8, range <<= 8;
		}
	}

	void End() override {
		for (int i = 0; i < 8; i++) {
            *data->pos++ = (uint8_t)(low >> 56);
			low <<= 8;
		}
	}
};

class CRangeDecoder : public CBasicRangeCoder {
public:
	CRangeDecoder() : buffer(0) { }

	void Start() override {
		buffer = 0;
		for (uint32_t i = 1; i <= 8; ++i) {
			buffer |= (uint64_t)*data->pos << (64 - i * 8);
            ++data->pos;
		}

		low = 0;
		range = Mask32;
	}

	inline uint32_t GetCumulativeFreq(uint32_t totalFreq) {
		assert(totalFreq != 0);
		return (uint32_t) (buffer / (range /= totalFreq));
	}

	void UpdateFrequency(uint32_t symFreq, uint32_t lowEnd, uint32_t /*totalFreq_*/) {
		uint32_t r = lowEnd * range;
		buffer -= r;
		low += r;
		range *= symFreq;

		while (range < TopValue) {
			if ( (uint8_t)((low ^ (low + range)) >> 56) )
				range = (((uint32_t)(low) | (TopValue - 1)) - (uint32_t)(low));

			buffer = (buffer << 8) + *data->pos++;
			low <<= 8, range <<= 8;
		}
	}

	void End() override {}

private:
	uint64_t buffer;
};


class CSimpleModel {
public: 
	CSimpleModel(uint32_t _step_size = 1) : n_symbols(0), stats(nullptr), step_size(_step_size)
	{
	};

	~CSimpleModel() {
		if (stats)
			delete[] stats;
	};

	CSimpleModel(const CSimpleModel &c) = delete;
	CSimpleModel& operator=(const CSimpleModel&) = delete;

	void Init(uint32_t _n_symbols, int* _init_stats, uint32_t _max_total, uint32_t _step_size) {
		step_size = _step_size;

		if (stats) {
			if (n_symbols != _n_symbols) {
				delete[] stats;
				n_symbols = _n_symbols;
				stats = new uint32_t[n_symbols];
			}
		} else {
			n_symbols = _n_symbols;
			stats = new uint32_t[n_symbols];
		}

		max_total = _max_total;

		if (_init_stats)
			for (uint32_t i = 0; i < n_symbols; ++i)
				stats[i] = _init_stats[i];
		else
			std::fill_n(stats, n_symbols, 1);

		total = std::accumulate(stats, stats+n_symbols, 0u);
		Rescale();
	}

	void Init(const CSimpleModel& c) {
		n_symbols = c.n_symbols;
		max_total = c.max_total;
		step_size = c.step_size;

		if (stats)
			delete[] stats;

		stats = new uint32_t[n_symbols];
		std::copy_n(c.stats, n_symbols, stats);
		total = std::accumulate(stats, stats + n_symbols, 0u);
	}

	void GetFreq(int symbol, int& sym_freq, int& left_freq, int& totf) {
		left_freq = 0;

		switch (symbol) {
			case 4: left_freq += stats[3];
			case 3: left_freq += stats[2];
			case 2: left_freq += stats[1];
			case 1: left_freq += stats[0];
			case 0: break;
			default:
				for (int i = 0; i < symbol; ++i)
					left_freq += stats[i];
		}

		sym_freq = stats[symbol];
		totf = total;
	}

	void Update(int symbol) {
		stats[symbol] += step_size;
		total += step_size;

		if (total >= max_total)
			Rescale();
	}

	int GetSym(int left_freq) const {
		int t = 0;

		for (uint32_t i = 0; i < n_symbols; ++i) {
			t += stats[i];
			if (t > left_freq)
				return i;
		}

		return -1;
	}

	inline uint32_t GetTotal() const { return total; }
	inline uint32_t* GetStats() const { return stats; }

	void SetStats(uint32_t* stats_to_set) {
		total = 0;
		for (uint32_t i = 0; i < n_symbols; ++i)
			total += stats[i] = stats_to_set[i];
	}

private:
    uint32_t n_symbols;
	uint32_t max_total;
	uint32_t* stats;
	uint32_t total;
	uint32_t step_size;

	void Rescale() {
		// std::cerr << "rescaling: " << total << ">=" << max_total << std::endl;
		while (total >= max_total) {
			total = 0;
			for (uint32_t i = 0; i < n_symbols; ++i) {
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
class CRangeCoderModel {
public:
	CRangeCoderModel(std::shared_ptr<CBasicRangeCoder> rcb, int _no_symbols, int _lg_totf, int _rescale, int* _init, uint32_t _step_size, bool _compress) :
		no_symbols(_no_symbols), lg_totf(_lg_totf), totf(1 << _lg_totf), rescale(_rescale), step_size(_step_size), compress(_compress)
	{
		simple_model.Init(no_symbols, _init, rescale, step_size);

		if (compress) {
			rce = std::static_pointer_cast<CRangeEncoder>(rcb);
			rcd = nullptr;
		} else {
			rce = nullptr;
			rcd = std::static_pointer_cast<CRangeDecoder>(rcb);
		}
	}

	CRangeCoderModel(const CRangeCoderModel &c) {
		simple_model.Init(c.simple_model);
		rce = c.rce;
		rcd = c.rcd;

		no_symbols = c.no_symbols;
		lg_totf = c.lg_totf;
		totf = c.totf;
		rescale = c.rescale;
		compress = c.compress;
		step_size = c.step_size;
	}

	~CRangeCoderModel() { }

/*	void Reset()
	{
		simple_model.Init(no_symbols, rescale);
	}*/

	void Encode(int x) {
		int syfreq, ltfreq;
		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	int Decode() {
		int syfreq, ltfreq;

		totf = simple_model.GetTotal();
		ltfreq = rcd->GetCumulativeFreq(totf);

		int x = simple_model.GetSym(ltfreq);

/*		if (x < 0)
		{
			cout << "Stats:\n";
			for (int i = 0; i < 4; ++i)
				cout << i << ": " << simple_model.stats[i] << "\n";
			cout << "totf: " << totf << endl;
			cout << "ltfreq: " << ltfreq << endl;
		}*/

		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rcd->UpdateFrequency(syfreq, ltfreq, totf);
		simple_model.Update(x);

		return x;
	}

	CSimpleModel* GetSimpleModel() { return &simple_model; }
	void Init(int *init) { simple_model.Init(no_symbols, init, rescale, step_size); }

public:
	std::shared_ptr<CRangeEncoder> rce;
	std::shared_ptr<CRangeDecoder> rcd;

	CSimpleModel simple_model;

	int no_symbols;
	int lg_totf;
	int totf;
	int rescale;
	uint32_t step_size;
	bool compress;
};

template<typename MODEL> 
class CContextHM {
public:
	typedef struct {
		uint64_t ctx;
		MODEL *rcm;
	} item_t;

	typedef uint64_t key_type;
	typedef MODEL* value_type;

private:
	double max_fill_factor;

	size_t size;
	size_t filled;
	item_t *data;
	size_t allocated;
	size_t size_when_restruct;
	size_t allocated_mask;

	size_t ht_memory;
	size_t ht_total;
	size_t ht_match;

	void restruct(void) {
		item_t *old_data = data;
		size_t old_allocated = allocated;

		allocated *= 2;
		size = 0;
		filled = 0;

		allocated_mask = allocated - 1ULL;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		ht_memory += allocated * sizeof(item_t);

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].rcm != nullptr)
				insert(old_data[i].ctx, old_data[i].rcm);

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);
	}

	inline size_t hash(uint64_t ctx) const { return (0x9e3779b97f4a7c13ULL * ctx) & allocated_mask; }

public:
	CContextHM() {
		ht_memory = 0;
		ht_total = 0;
		ht_match = 0;

		allocated = 16;
		allocated_mask = allocated - 1;

		size = 0;
		filled = 0;
		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		max_fill_factor = 0.4;

		ht_memory += allocated * sizeof(item_t);

		size_when_restruct = (size_t)(allocated * max_fill_factor);
	}

	~CContextHM() {
		if (data == nullptr)
			return;

		for (size_t i = 0; i < allocated; ++i) {
			if (data[i].rcm)
				delete data[i].rcm;
		}
		delete[] data;
	}

	inline size_t get_bytes() const { return ht_memory; }

	bool insert(uint64_t ctx, MODEL* rcm) {
		if (size >= size_when_restruct)
			restruct();

		size_t h = hash(ctx);

		if (data[h].rcm != nullptr) {
			do {
				h = (h + 1) & allocated_mask;
			} while (data[h].rcm != nullptr);
		}

		if (data[h].rcm == nullptr)
			++size;

		++filled;

		data[h].ctx = ctx;
		data[h].rcm = rcm;

		return true;
	}

	MODEL* find(uint64_t ctx) {
		size_t h = hash(ctx);

		if (data[h].rcm == nullptr)
			return nullptr;

		if (data[h].ctx == ctx)
			return data[h].rcm;

		h = (h + 1) & allocated_mask;

		while (data[h].rcm != nullptr) {
			if (data[h].ctx == ctx)
				return data[h].rcm;
			h = (h + 1) & allocated_mask;
		}

		return nullptr;
	}

	void Prefetch(uint64_t ctx) {
		size_t h = hash(ctx);

#ifdef _WIN32
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
#endif
	}

	size_t get_size(void) const { return filled; }
}; 

class djinn_gt_ctx {
public:
	djinn_gt_ctx() : ctx_prefix(0), ctx_symbol(0) {}

	inline CRangeCoderModel* GetEncodeModel(uint64_t ctx, uint32_t no_symbols, uint32_t max_log_counter) {
        CRangeCoderModel* p = rce_coders->find(ctx);

        if (p == nullptr) {
            rce_coders->insert(ctx, p = new CRangeCoderModel(rce, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, 1, true));
            p->rce->data = data;
        }

        return p;
    }

	inline CRangeCoderModel* GetEncodeModel(uint64_t ctx, uint32_t no_symbols, uint32_t max_log_counter, uint32_t add) {
        CRangeCoderModel* p = rce_coders->find(ctx);

        if (p == nullptr) {
            rce_coders->insert(ctx, p = new CRangeCoderModel(rce, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, add, true));
            p->rce->data = data;
        }

        return p;
    }

    inline CRangeCoderModel* GetDecodeModel(uint64_t ctx, uint32_t no_symbols, uint32_t max_log_counter) {
        CRangeCoderModel* p = rcd_coders->find(ctx);

        if (p == nullptr) {
            rcd_coders->insert(ctx, p = new CRangeCoderModel(rcd, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, 1, false));
            p->rcd->data = data;
        }

        return p;
    }

	inline CRangeCoderModel* GetDecodeModel(uint64_t ctx, uint32_t no_symbols, uint32_t max_log_counter, uint32_t add) {
        CRangeCoderModel* p = rcd_coders->find(ctx);

        if (p == nullptr) {
            rcd_coders->insert(ctx, p = new CRangeCoderModel(rcd, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, add, false));
            p->rcd->data = data;
        }

        return p;
    }

    int EncodeRunLength(uint8_t symbol, uint32_t len) {
        // Encode symbol
        CRangeCoderModel* rc_sym = GetEncodeModel(ctx_symbol + ctx_symbol_flag, 2, 15);
        rc_sym->Encode(symbol);
        ctx_symbol <<= 1;
        ctx_symbol += (symbol & 1);
        ctx_symbol &= ctx_symbol_mask;

        rce_coders->Prefetch(ctx_symbol + ctx_symbol_flag);

        ctx_prefix <<= 1;
        ctx_prefix += (uint64_t)(symbol & 1);
        ctx_prefix &= ctx_prefix_mask;

        // Encode run length
        CRangeCoderModel* rc_prefix = GetEncodeModel(ctx_prefix + ctx_prefix_flag, 11, 11);

        uint32_t prefix = ilog22(len);

        ctx_prefix <<= 4;
        ctx_prefix += (uint64_t)prefix;
        ctx_prefix &= ctx_prefix_mask;

        rce_coders->Prefetch(ctx_prefix + ctx_prefix_flag);

        if (prefix < 2) {
            rc_prefix->Encode(prefix);
        } 
		else if (prefix < 9) {
			assert(len < 512);
            rc_prefix->Encode(prefix);
            uint64_t ctx_suf = ctx_suffix_flag;
            ctx_suf += ((uint64_t)symbol); // 1 bit symbol + 4 bit log2(run)
            ctx_suf += (uint64_t) prefix << 1;
            uint32_t max_value_for_this_prefix = 1u << (prefix - 1);
            CRangeCoderModel* rc_suffix = GetEncodeModel(ctx_suf, max_value_for_this_prefix, 11);
            rc_suffix->Encode(len - max_value_for_this_prefix);
        } 
		else if (prefix < 16) {
			assert(len < 65536);
            rc_prefix->Encode(9); // flag for medium value

            uint64_t ctx_medium1 = ctx_medium_value1_flag;
            ctx_medium1 += (uint64_t)symbol;
            CRangeCoderModel* rc_l1 = GetEncodeModel(ctx_medium1, 256, 11);
            uint32_t lv1 = (len >> 8) & 255;
            rc_l1->Encode(lv1);

            uint64_t ctx_medium2 = ctx_medium_value2_flag;
            ctx_medium2 += (uint64_t)symbol;
            ctx_medium2 += (uint64_t)lv1 << 1;
            CRangeCoderModel* rc_l2 = GetEncodeModel(ctx_medium2, 256, 11);
            uint32_t lv2 = len & 255;
            rc_l2->Encode(lv2);
        } 
		else {
            rc_prefix->Encode(10); // flag for large value

            uint64_t ctx_large1 = ctx_large_value1_flag;
            ctx_large1 += ((uint64_t)symbol) << 16;
            CRangeCoderModel* rc_l1 = GetEncodeModel(ctx_large1, 256, 15);
            uint32_t lv1 = (len >> 16) & 0xff;
            rc_l1->Encode(lv1);

            uint64_t ctx_large2 = ctx_large_value2_flag;
            ctx_large2 += ((uint64_t)symbol) << 16;
            ctx_large2 += (uint64_t) lv1;
            CRangeCoderModel* rc_l2 = GetEncodeModel(ctx_large2, 256, 15);
            uint32_t lv2 = (len >> 8) & 0xff;
            rc_l2->Encode(lv2);

            uint64_t ctx_large3 = ctx_large_value3_flag;
            ctx_large3 += ((uint64_t)symbol) << 16;
            ctx_large3 += ((uint64_t) lv1) << 8;
            ctx_large3 += (uint64_t) lv2;
            CRangeCoderModel* rc_l3 = GetEncodeModel(ctx_large3, 256, 15);
            uint32_t lv3 = len & 0xff;
            rc_l3->Encode(lv3);
        }
    }

	int DecodeRunLength(uint8_t& symbol, uint32_t& len) {
		// Decode symbol
		auto rc_sym = GetDecodeModel(ctx_symbol + ctx_symbol_flag, 2, 15);
		symbol = (uint8_t) rc_sym->Decode();
		ctx_symbol <<= 1;
		ctx_symbol += (uint64_t) symbol;
		ctx_symbol &= ctx_symbol_mask;

		rcd_coders->Prefetch(ctx_symbol + ctx_symbol_flag);

		ctx_prefix <<= 1;
		ctx_prefix += (uint64_t)symbol;
		ctx_prefix &= ctx_prefix_mask;

		// Decode run length
		auto rc_p = GetDecodeModel(ctx_prefix + ctx_prefix_flag, 11, 9);

		uint32_t prefix = rc_p->Decode();

		if (prefix < 2) len = prefix;
		else if (prefix < 9) {
			uint64_t ctx_suf = ctx_suffix_flag;
			ctx_suf += ((uint64_t)symbol);
			ctx_suf += (uint64_t)prefix << 1;
			uint32_t max_value_for_this_prefix = 1u << (prefix - 1);

			auto rc_s = GetDecodeModel(ctx_suf, max_value_for_this_prefix, 9);
			len = max_value_for_this_prefix + rc_s->Decode();
			assert(len < 512);
		}
		else if (prefix == 9) {
			uint64_t ctx_large1 = ctx_medium_value1_flag;
			ctx_large1 += ((uint64_t)symbol);
			auto rc_l1 = GetDecodeModel(ctx_large1, 256, 10);
			uint32_t lv1 = rc_l1->Decode();

			uint64_t ctx_large2 = ctx_medium_value2_flag;
			ctx_large2 += ((uint64_t)symbol) ;
			ctx_large2 += (uint64_t) lv1 << 1;
			auto rc_l2 = GetDecodeModel(ctx_large2, 256, 10);
			uint32_t lv2 = rc_l2->Decode();

			len = (lv1 << 8) + lv2;
			assert(len < 65536);
			prefix = ilog22(len);
        } else {
			assert(prefix == 10);
			uint64_t ctx_large1 = ctx_large_value1_flag;
			ctx_large1 += ((uint64_t)symbol) << 16;
			auto rc_l1 = GetDecodeModel(ctx_large1, 256, 15);
			uint32_t lv1 = rc_l1->Decode();

			uint64_t ctx_large2 = ctx_large_value2_flag;
			ctx_large2 += ((uint64_t)symbol) << 16;
			ctx_large2 += (uint64_t) lv1;
			auto rc_l2 = GetDecodeModel(ctx_large2, 256, 15);
			uint32_t lv2 = rc_l2->Decode();

			uint64_t ctx_large3 = ctx_large_value3_flag;
			ctx_large3 += ((uint64_t)symbol) << 16;
			ctx_large3 += ((uint64_t) lv1) << 8;
			ctx_large3 += (uint64_t) lv2;
			auto rc_l3 = GetDecodeModel(ctx_large3, 256, 15);
			uint32_t lv3 = rc_l3->Decode();

			len = (lv1 << 16) + (lv2 << 8) + lv3;

			prefix = ilog22(len);
		}

		ctx_prefix <<= 4;
		ctx_prefix += (uint64_t)prefix;
		ctx_prefix &= ctx_prefix_mask;

		rcd_coders->Prefetch(ctx_prefix + ctx_prefix_flag);
	}

public:
	static constexpr uint64_t ctx_symbol_mask = (1 << 16) - 1;
	static constexpr uint64_t ctx_prefix_mask = (1 << 15) - 1;

    static constexpr uint64_t ctx_symbol_flag = 1ULL << 60;
	static constexpr uint64_t ctx_prefix_flag = 2ULL << 60;
	static constexpr uint64_t ctx_suffix_flag = 3ULL << 60;
	static constexpr uint64_t ctx_large_value1_flag = 4ULL << 60;
	static constexpr uint64_t ctx_large_value2_flag = 5ULL << 60;
	static constexpr uint64_t ctx_large_value3_flag = 6ULL << 60;
    static constexpr uint64_t ctx_medium_value1_flag = 7ULL << 60;
	static constexpr uint64_t ctx_medium_value2_flag = 8ULL << 60;
	
	uint64_t ctx_prefix;
	uint64_t ctx_symbol;
	
	std::shared_ptr<CContextHM<CRangeCoderModel>> rce_coders; // table of encode models
 	std::shared_ptr<CContextHM<CRangeCoderModel>> rcd_coders; // table of decode models

	std::shared_ptr<CRangeEncoder> rce; // shared RC encoder
    std::shared_ptr<CRangeDecoder> rcd; // shared RC decoder

    std::shared_ptr<djinn_ctx_buf_t> data;
};