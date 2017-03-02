#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <list>
using namespace std;

#ifndef PI2 /* PI/2 */
#define PI2 (2*atan((float)1))
#endif

/*directions
 * 3|2|1
 * 4|0|8
 * 5|6|7 */
extern int nextr[8];
extern int nextc[8];

typedef enum
{
	ZERO, /* zero cats do not accept zero category */
	FL, /* flat */
	PK, /* peak, summit */
	RI, /* ridge */
	SH, /* shoulder */
	CV, /* convex slope */
	SL, /* slope */
	CN, /* concave slope */
	FS, /* footslope */
	VL, /* valley */
	PT, /* pit, depression */
	__, /* error */
	CNT /* counter */
} FORMS;

typedef enum
{
	ZERON=0,
	FLN=1,//FL
	PKN=2,//PK
	RIN=3,//RI
	SHN=4,//SH
	CVN=5,//CV
	SLN=6,//SL
	CNN=7,//CN
	FSN=8,//FS
	VLN=9,//VL
	PTN=10,//PT
	__N=100,//__
	CNTN=11//CNT
} FORMS_NUM;

/*typedef struct 
{
	int i; //�ĸ�����
	int col; //x
	int row;//y
	float elevation;//��ʱ�洢��
	float height;//��ʱ�洢�߲�
	float cur_distance; //�����������ĵ����//  [1/17/2015 k]
} FeaturePoint;*/

typedef struct 
{
	int i;
	int col;
	int row;
	double angel;
	float elevation;
	float height;
	float distance;

}GeoPoint;

typedef struct 
{
	int num_positives;
	int num_negatives;
	unsigned char positives;
	unsigned char negatives;
	int pattern[8];// ʵ����ÿ������� ֵ 1 -1 0 Ӧ��ֻ��8��ֵ
	//float elevation[8];//ͬ��
	//double distance[8];//ͬ��
	//double x[8],y[8]; /* cartesian coordinates of geomorphon */
	// [11/20/2014 k] add by kangx
	GeoPoint feature[8]; //8�������feature 
} Pattern;

int extern determine_form(int num_minus, int num_plus);
int extern determine_form_num(int num_minus, int num_plus);
double extern caldist(int row1,int col1,int row2,int col2,double cell_res);
//void extern DoDPCalcu(vector <GeoPoint> Points, double Tolerance,vector <GeoPoint> &output);
int extern DPCalcuOnetime(vector<GeoPoint> Points, int firstpoint, int lastpoint, double Tolerance, vector<int> &PointToKeep);
void ParaCalcu(vector <GeoPoint> Points,int firstpoint, int lastpoint,double* k,double* b);
double DisCalcu(vector <GeoPoint> Points, int thispoint, double k,double b);	