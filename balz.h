// This program is extracted from BALZ 1.15 of Ilya Muraviev. 
// Source: http://balz.sourceforge.net/.
// ROLZ is removed and adaptive processing of bits is left.
// Each bit is predicted based on one preceding byte and two
// preceding bytes. One of these prediction is passed to encoder
// and decoder. 
//
//
// #define _CRT_DISABLE_PERFCRIT_LOCKS // for vc8 and later
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <memory.h>
#include <cstdint>

#define IDENT 0xba
// uint8_t* g_buffer = 0;
// int current_byte = 0;

///////auxiliary functions///////////////////////////////
// int getFileSize(char* fileName) {
// 	FILE* f;
// 	if ((f = fopen(fileName, "rb")) == NULL)
// 		return -1;
// 	fseek(f, 0, SEEK_END);
// 	int bytes = ftell(f);
// 	fclose(f);
// 	return bytes;
// }
//end

//Input/Output functions
// void writeByte(uint8_t byte) {
// 	g_buffer[current_byte++] = byte;
// }

// uint8_t readByte() {
// 	uint8_t byte = g_buffer[current_byte];
// 	++current_byte;
// 	return byte;
// }
//end

//This is predictor of Ilya Muraviev
class TPredictor {
private:
	uint16_t p1, p2;
public:
	TPredictor(): p1(1 << 15), p2(1 << 15) {} 
	~TPredictor() {}
	int P()	{
		return (p1 + p2); 
	}
	void Update(int bit) { 
		if (bit) {
			p1 += uint16_t(~p1) >> 3; 
			p2 += uint16_t(~p2) >> 6; 
		}
		else {
			p1 -= p1 >> 3; 
			p2 -= p2 >> 6; 
		}
	}

    void reset() {
        p1 = 1 << 15;
        p2 = 1 << 15;
    }
};

//arithmetic encoder and decoder derived by Ilya Muraviev from matt mahoney's fpaq0
class TBase {
public:
	TBase() : buffer(nullptr), dat(nullptr) {}
	uint8_t* buffer;
    uint8_t* dat;
};

class TEncoder : public TBase {
private:
	uint32_t x1, x2;
public:
	TEncoder(): x1(0), x2(-1) {} 

	void Encode(int probability, int bit) { 
		const uint32_t xmid = x1 + (uint32_t)((uint64_t(x2 - x1) * uint64_t (probability)) >> 17);
		if (bit) x2 = xmid;
		else x1 = xmid + 1;

		while ((x1 ^ x2) < (1 << 24)) {
			*dat = (uint8_t)(x2 >> 24);
            ++dat;
			x1 <<= 8;
			x2 = (x2 << 8) + 255;
		}
	}

	void Flush() { 
		for (int i=0; i<4; i++) {
			// writeByte((uint8_t)(x2 >> 24));
            *dat = (uint8_t)(x2 >> 24);
            ++dat;
			x2 <<= 8;
		}
	}

    void reset() {
        x1 = 0;
        x2 = -1;
        dat = buffer;
    }
};

class TDecoder : public TBase {
private:
	uint32_t x1, x2, x;
public:
	TDecoder(): x1(0), x2(-1) {} 

	void Init() { 
		for (int i=0; i<4; i++) {
			x = (x << 8) + *dat;
			++dat;
		}
	}

	int Decode(int P) {  

		const uint32_t xmid = x1 + (uint32_t)((uint64_t(x2 - x1) * uint64_t(P)) >> 17);
		int bit = (x <= xmid);
		if (bit) x2 = xmid;
		else x1 = xmid + 1;

		while ((x1 ^ x2) < (1 << 24)) { 
			x1 <<= 8;
			x2 = (x2 << 8) + 255;
			x = (x << 8) + *dat;
			++dat;
		}
		return (bit);
	}
};

class TPPM {
public:
	TPredictor** p1; //one byte context
	TPredictor** p2; //two bytes context
	uint32_t m_context;

public:

	TPPM() {
		p1 = 0;
		p1 = (TPredictor**)malloc(256 * sizeof(*p1));
		for (int i=0; i<256; ++i) p1[i] = new TPredictor[256]; 
		
		p2 = 0;
		p2 = (TPredictor**)malloc(65536 * sizeof(*p2));
		for (int i=0; i<65536; ++i) p2[i] = new TPredictor[256]; 

		m_context = 0;
	}
	~TPPM() {
		if (p1) {
			for (int i=0; i<256; ++i) {
				delete[] p1[i];
			}
			free(p1);
		}

		if (p2) {
			for (int i=0; i<65536; ++i) {
				delete[] p2[i];
			}
			free(p2);
		}

	}

	TEncoder encoder;
	TDecoder decoder;

    void reset() {
        m_context = 0;
        encoder.reset();
        for (int i=0; i<256; ++i){
			for (int j=0; j<256; ++j) {
				p1[i][j].reset();
			}
		}
        for (int i=0; i<65536; ++i){
			for (int j=0; j<256; ++j) {
				p2[i][j].reset();
			}
		}
    }

	void Encode(int value) { 
		for (int i=7, j=1; i>=0; i--) {  
			int bit = (value >> i) & 1;

			int probability1 = p1[m_context & 0xff][j].P();
			int probability = probability1;

			int probability2 = p2[m_context & 0xffff][j].P();
			if (abs(probability2 - 0xffff) > abs(probability1 - 0xffff)) {
				// std::cerr << "choosing higher" << std::endl;
				probability = probability2;
			} 
			// else {
			// 	std::cerr << "choosing lower" << std::endl;
			// }

			encoder.Encode(probability, bit);
			p1[m_context & 0xff][j].Update(bit);
			p2[m_context & 0xffff][j].Update(bit); 

			j += j + bit;
		}

		m_context <<= 8;
		m_context |= value;
	}

	int Decode() { 
		
		int value = 1;
		
		do {
			int probability1 = p1[m_context & 0xff][value].P();
			int probability = probability1;
			int probability2 = p2[m_context & 0xffff][value].P();
			if (abs(probability2 - 0xffff) > abs(probability1 - 0xffff)) {
				probability = probability2;
			}
			int bit = decoder.Decode(probability);
			p1[m_context & 0xff][value].Update(bit);

			p2[m_context & 0xffff][value].Update(bit); 
			value += value + bit;
		} while (value < 256); 

		int ret = value - 256;
		m_context <<= 8;
		m_context |= ret;


		return (ret);
	}
};


// encode in to out
// compressed file format:
// 1-byte - identification byte (0xba)
// 8-bytes - uncompressed size
// ?-bytes - compressed data
// void encode(uint8_t* data, int data_size) {
// 	TPPM ppm;
//     ppm.encoder.buffer = new uint8_t[10000000];
//     ppm.encoder.dat = ppm.encoder.buffer;

// 	uint64_t size=0; // uncompressed size
// 	// for (int i = 0; i < (1 + sizeof(size)); i++) {
// 	// 	writeByte(0); //we reserve space and write size when finished
// 	// }

// 	// exetransform(1, data, data_size); // perform a special exe transformation
// 	int i=0;
// 	while ((i < 2) && (i < data_size)) {
// 		ppm.Encode(data[i++], 0);
// 	}

// 	while (i < data_size) {
// 		ppm.Encode(data[i], data[i-1]);
// 		i++;
// 	}
// 	size += data_size; 

// 	ppm.encoder.Flush(); 
	
//     // //update signature and size
// 	// g_buffer[0] = IDENT;
// 	// for (int i=1; i<=sizeof(size); ++i) {
// 	// 	g_buffer[i] = (uint8_t)((size >> (8 * (i-1))) & 0xff);
// 	// }
// }

// void decode(uint8_t* result, int& result_size) {
// 	TPPM ppm;
// 	// check identification byte
// 	current_byte = 0;
// 	uint8_t byte = readByte();
// 	if (byte != IDENT) {
// 		printf("Bad file format\n");
// 		exit(1);
// 	}
// 	//read size of uncompressed data
// 	uint64_t size = 0x00;
// 	for (int i=0; i<sizeof(size); ++i) {
// 		size |= (readByte() << (8 * i));
// 	}
// 	if (size < 0) {
// 		fprintf(stderr, "size error\n");
// 		exit(1);
// 	}
// 	ppm.decoder.Init(); 

// 	result_size = 0;
// 	while ((result_size < 2) && (result_size < size)) { 
// 		const int value = ppm.Decode(0);
// 		result[result_size++] = value;
// 	}

// 	while (result_size < size) {
// 		const int value = ppm.Decode(result[result_size-1]); 
// 		result[result_size++] = value;
// 	}
// 	exetransform(0, result, result_size); 
// }

// int main() {

// 	printf("Extraction from balz v1.15 by ilia muraviev with removed ROLZ.\n");
// 	//
// 	char fileName[] = "/Users/Mivagallery/Downloads/1kgp3_chr20_subset.ubcf";
// 	int data_size = getFileSize(fileName);
// 	if (data_size < 0) {
// 		printf("File not found\n");
// 		return -1;
// 	}

// 	uint8_t* data = (uint8_t*)malloc(data_size * sizeof(*data));
// 	FILE* f = fopen(fileName, "rb");
// 	fread(data, 1, data_size, f);
// 	fclose(f);

// 	//we need copy for testing because original data are changed
// 	uint8_t* copy = (uint8_t*)malloc(data_size);
// 	memcpy(copy, data, data_size);

// 	//Encoding
// 	clock_t start = clock();
// 	g_buffer = (uint8_t*)malloc(data_size + data_size/2);
// 	encode(data, data_size);
// 	int compressed_size = current_byte;
// 	clock_t end = clock();
// 	printf("Encoding done, time %2.3f s.\n", (double)(end - start)/CLOCKS_PER_SEC);
// 	//End

// 	//Decoding
// 	clock_t start1 = clock();
// 	int test_size = data_size + data_size/2;
// 	uint8_t* test = (uint8_t*)malloc(test_size);
// 	decode(test, test_size);
// 	clock_t end1 = clock();
// 	printf("Decoding done, time %2.3f s.\n", (double)(end1 - start1)/CLOCKS_PER_SEC);
// 	//End

// 	bool isOK = true;
// 	for (int i=0; i<data_size; ++i) {
// 		if (test[i] != copy[i]) {
// 			printf("Data mismatch %8d %8d %8d\n", i, test[i], copy[i]);
// 			isOK = false;
// 			break;
// 		}
// 	}
// 	if (isOK) {
// 		printf("Round trip is OK\n");
// 	}
// 	printf("The compressed size relative to original %f\n", double(compressed_size)/double(data_size));

// 	if (test)     free(test);
// 	if (g_buffer) free(g_buffer);
// 	if (copy)     free(copy);
// 	if (data)     free(data);

// 	return (0);
// }