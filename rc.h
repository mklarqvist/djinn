#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <memory>//shared_ptr

// debug
#include <iostream>

template<typename T>
uint32_t ilog22(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

struct djinn_ctx_buf_t {
    djinn_ctx_buf_t() : buffer(nullptr), pos(nullptr), capac(0) {}
    djinn_ctx_buf_t(size_t cap) : buffer(new uint8_t[cap]), pos(buffer), capac(cap) {}
    ~djinn_ctx_buf_t() { delete[] buffer; buffer = nullptr; }

    uint8_t* buffer;
    uint8_t* pos;
    size_t capac;
};

class CBasicRangeCoder
{
public:
    CBasicRangeCoder() : low(0), range(0){}

public:
	static const uint32_t TopValue = 0x00ffffffULL;
	static const uint64_t Mask64 = 0xff00000000000000ULL;
	static const uint32_t Mask32 = 0xffffffffULL;

    uint64_t low;
	uint32_t range;

    std::shared_ptr<djinn_ctx_buf_t> data;
};

class CRangeEncoder : public CBasicRangeCoder {
public:
	CRangeEncoder(){}

	void Start() {
		low = 0;
		range = Mask32;
	}

	void EncodeFrequency(uint32_t symFreq_, uint32_t cumFreq_, uint32_t totalFreqSum_) {
		assert(range > totalFreqSum_);
		range /= totalFreqSum_;
		low += range * cumFreq_;
		range *= symFreq_;

		while (range <= TopValue) {
			assert(range != 0);
			if ((low ^ (low + range)) & Mask64) {
				uint32_t r = (uint32_t)low;
				range = (r | TopValue) - r;
			}
            *data->pos = (uint8_t) (low >> 56);
            ++data->pos;
			low <<= 8, range <<= 8;
		}
	}

	void End() {
		for (int i = 0; i < 8; i++) {
            *data->pos = (uint8_t) (low >> 56);
            ++data->pos;
			low <<= 8;
		}
	}
};

class CRangeDecoder : public CBasicRangeCoder {
public:
	CRangeDecoder() : buffer(0) { }

	void Start()
	{
		buffer = 0;
		for (uint32_t i = 1; i <= 8; ++i)
		{
			buffer |= (uint64_t)*data->pos << (64 - i * 8);
            ++data->pos;
		}

		low = 0;
		range = Mask32;
	}

	uint32_t GetCumulativeFreq(uint32_t totalFreq_)
	{
		assert(totalFreq_ != 0);
		return (uint32_t) (buffer / (range /= totalFreq_));
	}

	void UpdateFrequency(uint32_t symFreq_, uint32_t lowEnd_, uint32_t /*totalFreq_*/)
	{
		uint32_t r = lowEnd_*range;
		buffer -= r;
		low += r;
		range *= symFreq_;

		while (range <= TopValue)
		{
			if ((low ^ (low + range)) & Mask64)
			{
				uint32_t r = (uint32_t)low;
				range = (r | TopValue) - r;
			}

			buffer = (buffer << 8) + *data->pos;
            ++data->pos;
			low <<= 8, range <<= 8;
		}
	}

	void End() {}

private:
	uint64_t buffer;
};


// *******************************************************************************************
//
// *******************************************************************************************
class CSimpleModel {
public: 
	CSimpleModel(uint32_t _adder = 1) : n_symbols(0), stats(nullptr), adder(_adder)
	{
	};

	~CSimpleModel()
	{
		if (stats)
			delete[] stats;
	};

	CSimpleModel(const CSimpleModel &c) = delete;
	CSimpleModel& operator=(const CSimpleModel&) = delete;

	void Init(uint32_t _n_symbols, int *_init_stats, uint32_t _max_total, uint32_t _adder)
	{
		adder = _adder;

		if (stats)
		{
			if (n_symbols != _n_symbols)
			{
				delete[] stats;
				n_symbols = _n_symbols;
				stats = new uint32_t[n_symbols];
			}
		}
		else
		{
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
		rescale();
	}

	void Init(const CSimpleModel &c)
	{
		n_symbols = c.n_symbols;
		max_total = c.max_total;
		adder = c.adder;

		if (stats)
			delete[] stats;

		stats = new uint32_t[n_symbols];
		std::copy_n(c.stats, n_symbols, stats);
		total = std::accumulate(stats, stats + n_symbols, 0u);
	}

	void GetFreq(int symbol, int &sym_freq, int &left_freq, int &totf)
	{
		left_freq = 0;

		switch (symbol)
		{
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

	void Update(int symbol)
	{
//		stats[symbol]++;
//		total++;
		stats[symbol] += adder;
		total += adder;

		if (total >= max_total)
			rescale();
	}

	int GetSym(int left_freq)
	{
		int t = 0;

		for (uint32_t i = 0; i < n_symbols; ++i)
		{
			t += stats[i];
			if (t > left_freq)
				return i;
		}

		return -1;
	}

	uint32_t GetTotal()
	{
		return total;
	}

	void Merge(uint32_t *stats_to_merge)
	{
		for (uint32_t i = 0; i < n_symbols; ++i)
		{
			stats[i] += stats_to_merge[i];
			total += stats_to_merge[i];
		}
	}

	void CompleteMerge()
	{
		rescale();
	}

	uint32_t *GetStats()
	{
		return stats;
	}

	void SetStats(uint32_t *stats_to_set)
	{
		total = 0;
		for (uint32_t i = 0; i < n_symbols; ++i)
			total += stats[i] = stats_to_set[i];
	}

private:
    uint32_t n_symbols;
	uint32_t max_total;
	uint32_t *stats;
	uint32_t total;
	uint32_t adder;

	void rescale()
	{
		while (total >= max_total)
		{
			total = 0;
			for (uint32_t i = 0; i < n_symbols; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}

//		if (adder > 1)
//			adder /= 2;
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
class CRangeCoderModel {
public:
	CRangeCoderModel(CBasicRangeCoder *rcb, int _no_symbols, int _lg_totf, int _rescale, int* _init, uint32_t _adder, bool _compress) :
		no_symbols(_no_symbols), lg_totf(_lg_totf), totf(1 << _lg_totf), rescale(_rescale), adder(_adder), compress(_compress)
	{
		simple_model.Init(no_symbols, _init, rescale, adder);

		if (compress)
		{
			rce = (CRangeEncoder*) (rcb);
			rcd = nullptr;
		}
		else
		{
			rce = nullptr;
			rcd = (CRangeDecoder*) (rcb);
		}
	}

	CRangeCoderModel(const CRangeCoderModel &c)
	{
		simple_model.Init(c.simple_model);
		rce = c.rce;
		rcd = c.rcd;

		no_symbols = c.no_symbols;
		lg_totf = c.lg_totf;
		totf = c.totf;
		rescale = c.rescale;
		compress = c.compress;
		adder = c.adder;
	}

	~CRangeCoderModel()
	{
	}

/*	void Reset()
	{
		simple_model.Init(no_symbols, rescale);
	}*/

	void Encode(int x)
	{
		int syfreq, ltfreq;
		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	int Decode()
	{
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

	CSimpleModel* GetSimpleModel()
	{
		return &simple_model;
	}

	void Init(int *init)
	{
		simple_model.Init(no_symbols, init, rescale, adder);
	}

public:
	CRangeEncoder* rce;
	CRangeDecoder* rcd;

	CSimpleModel simple_model;

	int no_symbols;
	int lg_totf;
	int totf;
	int rescale;
	uint32_t adder;
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

		allocated_mask = allocated - 1ull;
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

	inline size_t hash(uint64_t ctx) { return (0x9e3779b97f4a7c13ull * ctx) & allocated_mask; }

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

		for (size_t i = 0; i < allocated; ++i)
			if (data[i].rcm)
				delete data[i].rcm;
		delete[] data;
	}

	size_t get_bytes() const {
		return ht_memory;
	}

	// Mozna to przyspieszyc tak, zebyinsert wykorzystywal wiedze o tym gdzie skonczyl szukac find
	bool insert(uint64_t ctx, MODEL *rcm) {
		if (size >= size_when_restruct)
			restruct();

		size_t h = hash(ctx);

		if (data[h].rcm != nullptr)
		{
			do
			{
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

		while (data[h].rcm != nullptr)
		{
			if (data[h].ctx == ctx)
				return data[h].rcm;
			h = (h + 1) & allocated_mask;
		}

		return nullptr;
	}

	void prefetch(uint64_t ctx)
	{
		size_t h = hash(ctx);

#ifdef _WIN32
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
#endif
	}

	size_t get_size(void) const
	{
		return filled;
	}
}; 

class djinn_gt_ctx {
public:
    CRangeEncoder* rce;
    CRangeDecoder* rcd;

    const uint64_t context_symbol_flag = 1ull << 60;
	const uint64_t context_symbol_mask = 0xffff;

	const uint64_t context_prefix_mask = 0xfffff;
	const uint64_t context_prefix_flag = 2ull << 60;
	const uint64_t context_suffix_flag = 3ull << 60;
	const uint64_t context_large_value1_flag = 4ull << 60;
	const uint64_t context_large_value2_flag = 5ull << 60;
	const uint64_t context_large_value3_flag = 6ull << 60;

    const uint64_t context_medium_value1_flag = 7ull << 60;
	const uint64_t context_medium_value2_flag = 8ull << 60;
	
	uint64_t ctx_prefix;
	uint64_t ctx_symbol;
	
	typedef CContextHM<CRangeCoderModel> ctx_map_e_t;
	typedef CContextHM<CRangeCoderModel> ctx_map_d_t;

	ctx_map_e_t rce_coders;
	ctx_map_d_t rcd_coders;

	inline CRangeCoderModel* find_rce_coder(uint64_t ctx, uint32_t no_symbols, uint32_t max_log_counter) {
        CRangeCoderModel* p = rce_coders.find(ctx);

        if (p == nullptr) {
            rce_coders.insert(ctx, p = new CRangeCoderModel(rce, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, 1, true));
            p->rce->data = data;
        }

        return p;
    }

    inline ctx_map_d_t::value_type find_rcd_coder(uint64_t ctx, uint32_t no_symbols, uint32_t max_log_counter) {
        auto p = rcd_coders.find(ctx);

        if (p == nullptr) {
            rcd_coders.insert(ctx, p = new CRangeCoderModel(rcd, no_symbols, max_log_counter, 1 << max_log_counter, nullptr, 1, false));
            p->rcd->data = data;
        }

        return p;
    }

    int encode_run_len(uint8_t symbol, uint32_t len) {
        // Encode symbol
        auto rc_sym = find_rce_coder(ctx_symbol + context_symbol_flag, 4, 15);
        rc_sym->Encode(symbol);
        ctx_symbol <<= 4;
        ctx_symbol += symbol;
        ctx_symbol &= context_symbol_mask;

        rce_coders.prefetch(ctx_symbol + context_symbol_flag);

        ctx_prefix <<= 4;
        ctx_prefix += (uint64_t)symbol;
        ctx_prefix &= context_prefix_mask;

        // Encode run length
        auto rc_p = find_rce_coder(ctx_prefix + context_prefix_flag, 12, 10);

        uint32_t prefix = ilog22(len);

        ctx_prefix <<= 4;
        ctx_prefix += (uint64_t)prefix;
        ctx_prefix &= context_prefix_mask;

        rce_coders.prefetch(ctx_prefix + context_prefix_flag);

        if (prefix < 2)
            rc_p->Encode(prefix);
        else if (prefix < 10)
        {
            rc_p->Encode(prefix);
            uint64_t ctx_suf = context_suffix_flag;
            ctx_suf += ((uint64_t)symbol) << 8;
            ctx_suf += (uint64_t) prefix;
            uint32_t max_value_for_this_prefix = 1u << (prefix - 1);

            // std::cerr << len - max_value_for_this_prefix << "/" << max_value_for_this_prefix << std::endl;
            auto rc_s = find_rce_coder(ctx_suf, max_value_for_this_prefix, 15);
            rc_s->Encode(len - max_value_for_this_prefix);
        }
        else if (prefix < 16) {
            rc_p->Encode(10);		// flag for medium value

            uint64_t ctx_medium1 = context_medium_value1_flag;
            ctx_medium1 += ((uint64_t)symbol) << 16;
            auto rc_l1 = find_rce_coder(ctx_medium1, 256, 15);
            uint32_t lv1 = (len >> 8) & 0xff;
            rc_l1->Encode(lv1);

            uint64_t ctx_medium2 = context_medium_value2_flag;
            ctx_medium2 += ((uint64_t)symbol) << 16;
            ctx_medium2 += (uint64_t) lv1;
            auto rc_l2 = find_rce_coder(ctx_medium2, 256, 15);
            uint32_t lv2 = len & 0xff;
            rc_l2->Encode(lv2);
        }
        else {
            // std::cerr << "here" << std::endl;
            rc_p->Encode(11);		// flag for large value

            uint64_t ctx_large1 = context_large_value1_flag;
            ctx_large1 += ((uint64_t)symbol) << 16;
            auto rc_l1 = find_rce_coder(ctx_large1, 256, 15);
            uint32_t lv1 = (len >> 16) & 0xff;
            rc_l1->Encode(lv1);

            uint64_t ctx_large2 = context_large_value2_flag;
            ctx_large2 += ((uint64_t)symbol) << 16;
            ctx_large2 += (uint64_t) lv1;
            auto rc_l2 = find_rce_coder(ctx_large2, 256, 15);
            uint32_t lv2 = (len >> 8) & 0xff;
            rc_l2->Encode(lv2);

            uint64_t ctx_large3 = context_large_value3_flag;
            ctx_large3 += ((uint64_t)symbol) << 16;
            ctx_large3 += ((uint64_t) lv1) << 8;
            ctx_large3 += (uint64_t) lv2;
            auto rc_l3 = find_rce_coder(ctx_large3, 256, 15);
            uint32_t lv3 = len & 0xff;
            rc_l3->Encode(lv3);
        }
    }

    std::shared_ptr<djinn_ctx_buf_t> data;
};