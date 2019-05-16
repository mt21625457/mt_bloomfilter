#include "MtBloomfilter.h"





mt::mtBloomfilter::mtBloomfilter()
{
}

mt::mtBloomfilter::~mtBloomfilter()
{
}

bool mt::mtBloomfilter::init(uint32_t dwsend, uint32_t n, double p)
{
	return this->InitBloomFilter(dwsend, n, p);
}

bool mt::mtBloomfilter::add(const void * key, int len)
{
	return this->BloomFilter_Add(key,len);
}

bool mt::mtBloomfilter::find(const void * key, int len)
{
	return this->BloomFilter_Check(key,len);
}

void mt::mtBloomfilter::close()
{
	this->FreeBloomFilter();
}

uint64_t mt::mtBloomfilter::MurmurHash2_64_64(const void * key, int len, unsigned int seed)
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len / 8);

	while (data != end)
	{
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch (len & 7)
	{
	case 7: h ^= uint64_t(data2[6]) << 48;
	case 6: h ^= uint64_t(data2[5]) << 40;
	case 5: h ^= uint64_t(data2[4]) << 32;
	case 4: h ^= uint64_t(data2[3]) << 24;
	case 3: h ^= uint64_t(data2[2]) << 16;
	case 2: h ^= uint64_t(data2[1]) << 8;
	case 1: h ^= uint64_t(data2[0]);
		h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}

uint64_t mt::mtBloomfilter::MurmurHash2_64_32(const void * key, int len, unsigned int seed)
{
	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	unsigned int h1 = seed ^ len;
	unsigned int h2 = 0;

	const unsigned int * data = (const unsigned int *)key;

	while (len >= 8)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;

		unsigned int k2 = *data++;
		k2 *= m; k2 ^= k2 >> r; k2 *= m;
		h2 *= m; h2 ^= k2;
		len -= 4;
	}

	if (len >= 4)
	{
		unsigned int k1 = *data++;
		k1 *= m; k1 ^= k1 >> r; k1 *= m;
		h1 *= m; h1 ^= k1;
		len -= 4;
	}

	switch (len)
	{
	case 3: h2 ^= ((unsigned char*)data)[2] << 16;
	case 2: h2 ^= ((unsigned char*)data)[1] << 8;
	case 1: h2 ^= ((unsigned char*)data)[0];
		h2 *= m;
	};

	h1 ^= h2 >> 18; h1 *= m;
	h2 ^= h1 >> 22; h2 *= m;
	h1 ^= h2 >> 17; h1 *= m;
	h2 ^= h1 >> 19; h2 *= m;

	uint64_t h = h1;

	h = (h << 32) | h2;

	return h;
}

void mt::mtBloomfilter::CalcBloomFilterParam(uint32_t n, double p, uint32_t * pm, uint32_t * pk)
{
	/**
	 *  n - Number of items in the filter
	 *  p - Probability of false positives, float between 0 and 1 or a number indicating 1-in-p
	 *  m - Number of bits in the filter
	 *  k - Number of hash functions
	 *
	 *  f = ln(2) × ln(1/2) × m / n = (0.6185) ^ (m/n)
	 *  m = -1 * ln(p) × n / 0.6185
	 *  k = ln(2) × m / n = 0.6931 * m / n
	**/

	uint32_t m, k;

	// 计算指定假阳概率下需要的比特数
	m = (uint32_t)ceil(-1 * log(p) * n / 0.6185);
	m = (m - m % 64) + 64;                  // 8字节对齐

	// 计算哈希函数个数
	k = (uint32_t)(0.6931 * m / n);
	k++;

	*pm = m;
	*pk = k;
	return;
}

inline bool mt::mtBloomfilter::InitBloomFilter(uint32_t dwSeed, uint32_t dwMaxItems, double dProbFalse)
{
	if (this->m_pstBloomfilter != nullptr) {
		return false;
	}

	// 初始化内存结构，并计算BloomFilter需要的空间
	this->m_pstBloomfilter = new BaseBloomFilter;
	this->m_pstBloomfilter->dwMaxItems	= dwMaxItems;
	this->m_pstBloomfilter->dProbFalse	= dProbFalse;
	this->m_pstBloomfilter->dwSeed		= dwSeed;

	// 计算 m, k
	this->CalcBloomFilterParam(this->m_pstBloomfilter->dwMaxItems, this->m_pstBloomfilter->dProbFalse,
		&(this->m_pstBloomfilter->dwFilterBits), &(this->m_pstBloomfilter->dwHashFuncs));

	// 分配BloomFilter的存储空间
	this->m_pstBloomfilter->dwFilterSize = this->m_pstBloomfilter->dwFilterBits / BYTE_BITS;
	this->m_pstBloomfilter->pstFilter = new uchar_t[this->m_pstBloomfilter->dwFilterSize];

	// 初始化BloomFilter的内存
	memset(this->m_pstBloomfilter->pstFilter, 0, this->m_pstBloomfilter->dwFilterSize);

	//哈希结果数组，每个哈希函数一个
	this->m_pstBloomfilter->pdwHashPos = new uint32_t[this->m_pstBloomfilter->dwHashFuncs * sizeof(uint32_t)];

	return true;
}

inline void mt::mtBloomfilter::FreeBloomFilter()
{
	if (this->m_pstBloomfilter != nullptr) {
		this->m_pstBloomfilter->dwCount = 0;

		if (this->m_pstBloomfilter->pstFilter != nullptr) {
			delete[] this->m_pstBloomfilter->pstFilter;
			this->m_pstBloomfilter->pstFilter = nullptr;
		}

		if ( this->m_pstBloomfilter->pdwHashPos != nullptr) {
			delete[] this->m_pstBloomfilter->pdwHashPos;
			this->m_pstBloomfilter->pdwHashPos = nullptr;
		}

		delete this->m_pstBloomfilter;
		this->m_pstBloomfilter = nullptr;
	}
}

inline void mt::mtBloomfilter::bloom_hash_64(const void * key, int len)
{
	uint32_t dwFilterBits = this->m_pstBloomfilter->dwFilterBits;
	uint64_t hash1 = MurmurHash2_64_64(key, len, this->m_pstBloomfilter->dwSeed);
	uint64_t hash2 = MurmurHash2_64_64(key, len, MIX_UINT64(hash1));

	for (uint32_t i = 0; i < this->m_pstBloomfilter->dwHashFuncs; ++i) {
		this->m_pstBloomfilter->pdwHashPos[i] = (hash1 + i * hash2) % dwFilterBits;
	}
}

inline void mt::mtBloomfilter::bloom_hash_32(const void * key, int len)
{
}

inline bool mt::mtBloomfilter::BloomFilter_Add(const void * key, int len)
{
	if (this->m_pstBloomfilter == nullptr || key == nullptr || len <= 0) {
		return false;
	}

	if (this->m_pstBloomfilter->dwCount >= this->m_pstBloomfilter->dwMaxItems) {
		return false;
	}

	this->bloom_hash_64(key, len);
	for (uint32_t i = 0; i < this->m_pstBloomfilter->dwHashFuncs; ++i) {
		SETBIT(this->m_pstBloomfilter, this->m_pstBloomfilter->pdwHashPos[i]);
	}

	this->m_pstBloomfilter++;

	return true;
}

inline bool mt::mtBloomfilter::BloomFilter_Check(const void * key, int len)
{
	if ( this->m_pstBloomfilter == nullptr ) {
		return false;
	}

	this->bloom_hash_64(key, len);
	for (uint32_t i = 0; i < this->m_pstBloomfilter->dwHashFuncs; ++i) {
		if ( GETBIT(this->m_pstBloomfilter, this->m_pstBloomfilter->pdwHashPos[i]) == 0 ) {
			return false;
		}		
	}

	return true;
}

