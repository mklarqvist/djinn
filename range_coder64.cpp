#include "range_coder64.h"

#include <iostream>

namespace djinn {

RangeCoder64::RangeCoder64() :
	Low(0),
	Range((uint64_t)-1)
{
}

RangeEncoder64::RangeEncoder64() :
	Flushed(false)
{
}

RangeEncoder64::~RangeEncoder64()
{
	// if (!Flushed) Flush();
}

void RangeEncoder64::EncodeRange(uint8_t*& pos, uint32_t SymbolLow, uint32_t SymbolHigh, uint32_t TotalRange)
{
	Low += SymbolLow * (Range /= TotalRange);
	Range *= SymbolHigh - SymbolLow;

	while ((Low ^ (Low+Range)) < Top || Range < Bottom && ((Range= -Low & (Bottom-1)),1))
	{
		// Output.WriteByte(Low >> 56);
        // std::cout << (Low >> 56);
        *pos = Low >> 56;
        ++pos;
        Range <<= 8;
        Low <<= 8;
	}
}

void RangeEncoder64::Flush(uint8_t*& pos)
{
	// if(!Flushed)
	// {
		for(uint32_t i = 0; i < 8; ++i) {
			// Output.WriteByte(Low>>56);
            // std::cout << (Low >> 56);
            *pos = Low >> 56;
            ++pos;
			Low <<= 8;
		}
		// Flushed=true;
	// }
}

RangeDecoder64::RangeDecoder64() :
	Code(0)
{
	for(uint32_t i = 0; i < 8; ++i) {
        // Todo:
		// Code = (Code << 8) | Input.ReadByte();
	}
}

uint32_t RangeDecoder64::GetCurrentCount(uint32_t TotalRange)
{
	return (Code-Low) / (Range /= TotalRange);
}

void RangeDecoder64::RemoveRange(uint32_t SymbolLow, uint32_t SymbolHigh, uint32_t /*TotalRange*/)
{
	Low += SymbolLow * Range;
	Range *= SymbolHigh - SymbolLow;

	while ((Low ^ Low+Range) < Top || Range < Bottom && ((Range= -Low & Bottom-1),1))
	{
        // Todo
		// Code= Code<<8 | Input.ReadByte(), Range<<=8, Low<<=8;
	}
}

}