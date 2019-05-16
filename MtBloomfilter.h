/******************************************************************************
*  @MyClass.cpp  : ( ����Cassandra�е�BloomFilterʵ�֣�Hashѡ��MurmurHash2)			
*    ͨ��˫��ɢ�й�ʽ����ɢ�к������ο���http://hur.st/bloomfilter
*    Hash(key, i) = (H1(key) + i * H2(key)) % m
*
******************************************************************************
*  @author	:	mt
*  @date	:	2019/02/13
*  @version	:	0.1
******************************************************************************
*  ģ���б� :
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

// ע�⣬Ҫ��Add/Check��������������ʹ�� -O2 �����ϵ��Ż��ȼ�
//#define FORCE_INLINE __attribute__((always_inline))

#define BYTE_BITS           (8)
#define MIX_UINT64(v)       ((uint32_t)((v>>32)^(v)))

/*����λͼ*/
#define SETBIT(filter, n)   (filter->pstFilter[n/BYTE_BITS] |= (1 << (n%BYTE_BITS)))
#define GETBIT(filter, n)   (filter->pstFilter[n/BYTE_BITS] & (1 << (n%BYTE_BITS)))

/*�ֽڶ���*/


typedef unsigned __int64  uint64_t;
typedef unsigned int	  uint32_t;
typedef unsigned char	  uchar_t;

namespace mt {

enum hash_type {
	MURMURHASH2_64_64 = 0,
	MURMURHASH2_64_32
};


// BloomFilter�ṹ����
struct BaseBloomFilter
{
	uint32_t	dwCount;								// Add()�ļ���������MAX_BLOOMFILTER_N�򷵻�ʧ��
	uint32_t	dwMaxItems;								// n - BloomFilter�����Ԫ�ظ��� (������)
	double		dProbFalse;								// p - �������� (���������������֮һ��0.00001)
	uint32_t	dwFilterBits;							// m = ceil((n * log(p)) / log(1.0 / (pow(2.0, log(2.0))))); - BloomFilter�ı�����
	uint32_t	dwHashFuncs;							// k = round(log(2.0) * m / n); - ��ϣ��������

	uint32_t	dwSeed;									// MurmurHash������ƫ����
	uint32_t	dwFilterSize;							// dwFilterBits / BYTE_BITS

	uchar_t	 *	pstFilter;								// BloomFilter�洢ָ�룬  
	uint32_t *	pdwHashPos;								// �洢�ϴ�hash�õ���K��bitλ������(��bloom_hash���)
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

private: /* ��ϣ���� */
	/*64λƽ̨��64λɢ��*/
	uint64_t MurmurHash2_64_64(const void * key, int len, unsigned int seed);

	/*32λƽ̨��64λɢ��*/
	uint64_t MurmurHash2_64_32(const void * key, int len, unsigned int seed);

private:
	/*����BloomFilter�Ĳ���m,k*/
	inline void CalcBloomFilterParam(uint32_t n, double p, uint32_t *pm, uint32_t *pk);

	/*����Ŀ�꾫�Ⱥ����ݸ�������ʼ��BloomFilter�ṹ */
	inline bool InitBloomFilter(uint32_t dwSeed, uint32_t dwMaxItems, double dProbFalse);

	/*�ͷ�BloomFilter*/
	inline void FreeBloomFilter();

	/*64λ˫��ɢ�з�װ*/
	inline void bloom_hash_64(const void * key, int len);

	/*32λ˫��ɢ�з�װ*/
	inline void bloom_hash_32(const void * key, int len);

	/* ��BloomFilter������һ��Ԫ�� */
	inline bool BloomFilter_Add(const void * key, int len);

	/*���һ��Ԫ���Ƿ���bloomfilter�� */
	inline bool BloomFilter_Check(const void * key, int len);
};


}  // namespaec mt

#endif 