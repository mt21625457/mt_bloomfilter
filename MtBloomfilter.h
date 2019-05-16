/******************************************************************************
*  @MyClass.cpp  : ( 仿照Cassandra中的BloomFilter实现，Hash选用MurmurHash2)			
*    通过双重散列公式生成散列函数，参考：http://hur.st/bloomfilter
*    Hash(key, i) = (H1(key) + i * H2(key)) % m
*
******************************************************************************
*  @author	:	mt
*  @date	:	2019/02/13
*  @version	:	0.1
******************************************************************************
*  模块列表 :
*
*  Class       : 
*  Description :
*
******************************************************************************
*  Change History
* ----------------------------------------------------------------------------
*    Date    :  Ver. : Author : Description
* ----------------------------------------------------------------------------
* 2019/02/13 :       : mt : Create
*            :       :        :
******************************************************************************/

#ifndef __MT_MTBLOOMFILTER_H__
#define __MT_MTBLOOMFILTER_H__

#pragma pack(8)

extern "C" {
	#include <stdio.h>
	#include <stdlib.h>
	#include <stdint.h>
	#include <string.h>
	#include <math.h>
}



#define __MT_BLOOMFILTER_VERSION__ "0.1"
//#define __MGAIC_CODE__          (0x01464C42)

// 注意，要让Add/Check函数内联，必须使用 -O2 或以上的优化等级
//#define FORCE_INLINE __attribute__((always_inline))

#define BYTE_BITS           (8)
#define MIX_UINT64(v)       ((uint32_t)((v>>32)^(v)))

/*设置位图*/
#define SETBIT(filter, n)   (filter->pstFilter[n/BYTE_BITS] |= (1 << (n%BYTE_BITS)))
#define GETBIT(filter, n)   (filter->pstFilter[n/BYTE_BITS] & (1 << (n%BYTE_BITS)))

/*字节对齐*/


typedef unsigned __int64  uint64_t;
typedef unsigned int	  uint32_t;
typedef unsigned char	  uchar_t;

namespace mt {

enum hash_type {
	MURMURHASH2_64_64 = 0,
	MURMURHASH2_64_32
};


// BloomFilter结构定义
struct BaseBloomFilter
{
	uint32_t	dwCount;								// Add()的计数，超过MAX_BLOOMFILTER_N则返回失败
	uint32_t	dwMaxItems;								// n - BloomFilter中最大元素个数 (输入量)
	double		dProbFalse;								// p - 假阳概率 (输入量，比如万分之一：0.00001)
	uint32_t	dwFilterBits;							// m = ceil((n * log(p)) / log(1.0 / (pow(2.0, log(2.0))))); - BloomFilter的比特数
	uint32_t	dwHashFuncs;							// k = round(log(2.0) * m / n); - 哈希函数个数

	uint32_t	dwSeed;									// MurmurHash的种子偏移量
	uint32_t	dwFilterSize;							// dwFilterBits / BYTE_BITS

	uchar_t	 *	pstFilter;								// BloomFilter存储指针，  
	uint32_t *	pdwHashPos;								// 存储上次hash得到的K个bit位置数组(由bloom_hash填充)
} ;

class mtBloomfilter
{
public:
	mtBloomfilter();
	mtBloomfilter(uint32_t dwsend, uint32_t n, double p) = delete;
	mtBloomfilter(const mtBloomfilter &bf) = delete;
	~mtBloomfilter();

public:
	bool init(uint32_t dwsend, uint32_t n, double p);
	bool add(const void * key, int len);
	bool find(const void * key, int len);
	void close();

private:
	BaseBloomFilter * m_pstBloomfilter = nullptr;

private: /* 哈希函数 */
	/*64位平台的64位散列*/
	uint64_t MurmurHash2_64_64(const void * key, int len, unsigned int seed);

	/*32位平台的64位散列*/
	uint64_t MurmurHash2_64_32(const void * key, int len, unsigned int seed);

private:
	/*计算BloomFilter的参数m,k*/
	inline void CalcBloomFilterParam(uint32_t n, double p, uint32_t *pm, uint32_t *pk);

	/*根据目标精度和数据个数，初始化BloomFilter结构 */
	inline bool InitBloomFilter(uint32_t dwSeed, uint32_t dwMaxItems, double dProbFalse);

	/*释放BloomFilter*/
	inline void FreeBloomFilter();

	/*64位双重散列封装*/
	inline void bloom_hash_64(const void * key, int len);

	/*32位双重散列封装*/
	inline void bloom_hash_32(const void * key, int len);

	/* 向BloomFilter中新增一个元素 */
	inline bool BloomFilter_Add(const void * key, int len);

	/*检查一个元素是否在bloomfilter中 */
	inline bool BloomFilter_Check(const void * key, int len);
};


}  // namespaec mt

#endif 