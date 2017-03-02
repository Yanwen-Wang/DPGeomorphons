#include"slopeOperator.h"
#include"math.h"
#include"geomorphons.h"
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>


/*directions
 * 3|2|1
 * 4|0|8
 * 5|6|7 */
//int nextr[8] = {-1, -1, -1, 0, 1, 1, 1, 0 };
//int nextc[8] = {1, 0, -1, -1, -1, 0, 1, 1 };

void SlopeOperator::demLayer(RasterLayer<double> &layerD) 
{
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
  cellSize = _pDEMLayer->_pMetaData->cellSize;
  noData = _pDEMLayer->_pMetaData->noData;
    //pWorkBR = &_pDEMLayer->_pMetaData->_localworkBR;
  Configure(_pDEMLayer, false);
}

void SlopeOperator::slopeLayer(RasterLayer<double> &layerD) 
{
  _pSlopeLayer = &layerD;
  Configure(_pSlopeLayer,false);
}

bool SlopeOperator::isTermination()
{
	num--;
	if(num > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//�������
double caldist(int row1,int col1,int row2,int col2,double cell_res)
{
	double dist;
	dist = cell_res * sqrt((double)((col1-col2)*(col1-col2) +(row1-row2)*(row1-row2)));
	return dist;
}
//�����߲���
void ParaCalcu(vector<GeoPoint> Points,int firstpoint, int lastpoint,double *k,double *b)
{
	//LineParameters thisLine = new LineParameters();
	*k = ((double)Points[lastpoint].elevation-(double)Points[firstpoint].elevation)/((double)Points[lastpoint].distance-(double)Points[firstpoint].distance);
	*b = (double)Points[firstpoint].elevation - *k*(double)Points[firstpoint].distance;
	//cout<<"step1.1.1.1"<<" "<<firstpoint<<" "<<lastpoint<<endl;
	//cout<<"step1.1.1"<<" "<<Points[lastpoint].elevation<<" "<<Points[lastpoint].distance<<" "<<Points[firstpoint].elevation<<" "<<Points[firstpoint].distance<<endl;
	//cout<<"step1.1"<<" "<<k<<" "<<b<<endl;
}
//����ÿһ����ֱ�ߵĸ߲�
double DisCalcu(vector<GeoPoint> Points, int thispoint, double k,double b)
{
	//GeoPoint temPoint = new GeoPoint(0,0,0);
	double tem = 0; double distance;			
	tem = Points[thispoint].distance * k + b;			
	distance = fabs(Points[thispoint].elevation - tem);

	return distance;
}
//����DPһ�μ���
int DPCalcuOnetime(vector<GeoPoint> Points, int firstpoint, int lastpoint, double Tolerance, vector<int> &PointToKeep)
{
	double maxDistance = 0;
	int FarthestPoint = 0;
	double k,b;

	//����ÿһ����ֱ�ߵĸ߲�
	ParaCalcu(Points, firstpoint, lastpoint,&k,&b);
	//cout<<"step1"<<k<<" "<<b<<endl;
	for (int thispoint = firstpoint;thispoint<lastpoint;thispoint++)
	{
		double distance = DisCalcu(Points,thispoint,k,b);
		//cout<<"distance="<<distance<<" k="<<k<<" b="<<b<<endl;
		//cout<<"step2"<<thispoint<<" "<<distance<<endl;
		if (distance > maxDistance)
		{
			maxDistance = distance;
			FarthestPoint = thispoint;
		}
	}
	//cout<<firstpoint<<" "<<lastpoint<<" "<<FarthestPoint<<endl;
	/*if (maxDistance > Tolerance && FarthestPoint != 0)
	{*/
		//��������㱣��
		PointToKeep.push_back(FarthestPoint);
		//�ݹ����
		//DPCalcu(Points, firstpoint, FarthestPoint, Tolerance,PointToKeep);
		//DPCalcu(Points, FarthestPoint, lastpoint, Tolerance, PointToKeep);
	//}
	//else
		//cout<<"FarthestPoint="<<FarthestPoint<<endl;
		return FarthestPoint;
}

//����patten����
int determine_form(int num_minus, int num_plus)
{
/* determine form according number of positives and negatives
 * simple approach to determine form pattern */
	const FORMS forms[9][9] = 
	{
/* minus ------------- plus ----------------*/
/*       0   1   2   3   4   5   6   7   8  */
/* 0 */ {FL, FL, FL, FS, FS, VL, VL, VL, PT},
/* 1 */ {FL, FL, FS, FS, FS, VL, VL, VL, __},
/* 2 */ {FL, SH, SL, SL, CN, CN, VL, __, __},
/* 3 */ {SH, SH, SL, SL, SL, CN, __, __, __},
/* 4 */ {SH, SH, CV, SL, SL, __, __, __, __},
/* 5 */ {RI, RI, CV, CV, __, __, __, __, __},
/* 6 */ {RI, RI, RI, __, __, __, __, __, __},
/* 7 */ {RI, RI, __, __, __, __, __, __, __},
/* 8 */ {PK, __, __, __, __, __, __, __, __},
    };

/* legend:
  FL,  flat
  PK,  peak, summit
  RI,  ridge
  SH,  shoulder
  CV,  convex
  SL,  slope
  CN,  concave
  FS,  footslope
  VL,  valley
  PT,  pit, depression
  __  error, impossible
*/
 return forms[num_minus][num_plus];
}

int determine_form_num(int num_minus, int num_plus)
{
/* determine form according number of positives and negatives
 * simple approach to determine form pattern */
	const FORMS_NUM forms_num[9][9] = 
	{
		/* minus ------------- plus ----------------*/
		/*       0   1   2   3   4   5   6   7   8  */
		/* 0 */ {FLN, FLN, FLN, FSN, FSN, VLN, VLN, VLN, PTN},
		/* 1 */ {FLN, FLN, FSN, FSN, FSN, VLN, VLN, VLN, __N},
		/* 2 */ {FLN, SHN, SLN, SLN, CNN, CNN, VLN, __N, __N},
		/* 3 */ {SHN, SHN, SLN, SLN, SLN, CNN, __N, __N, __N},
		/* 4 */ {SHN, SHN, CVN, SLN, SLN, __N, __N, __N, __N},
		/* 5 */ {RIN, RIN, CVN, CVN, __N, __N, __N, __N, __N},
		/* 6 */ {RIN, RIN, RIN, __N, __N, __N, __N, __N, __N},
		/* 7 */ {RIN, RIN, __N, __N, __N, __N, __N, __N, __N},
		/* 8 */ {PKN, __N, __N, __N, __N, __N, __N, __N, __N},
	};

/* legend:
  FL,  flat
  PK,  peak, summit
  RI,  ridge
  SH,  shoulder
  CV,  convex
  SL,  slope
  CN,  concave
  FS,  footslope
  VL,  valley
  PT,  pit, depression
  __  error, impossible
*/
 return forms_num[num_minus][num_plus];
}

bool SlopeOperator::Operator(const CellCoord &coord)
{
	CellSpace<double> &dem = *(_pDEMLayer->cellSpace());//����ͼ���դ������
	
	CellSpace<double> &slope = *(_pSlopeLayer->cellSpace());//���ͼ���դ������
	  
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);//���������ļ�

	int nextr[8] = {-1, -1, -1, 0, 0, 1, 1, 1 };
	int nextc[8] = {-1, 0, 1, -1, 1, -1, 0, 1 };
	
	//coord�Ƿ������ڵ�����ԭ�㣬Ҳ������Ҫ�����¶ȵĵ�
	int cur_Row = coord.iRow();//��ǰ��row
	int cur_Col = coord.iCol();//��ǰ��cols

	//neighborhood�ļ�һ����ʮ���㣬ǰ�Ÿ����ǲ����㣬��ʮ���Ǵ��ڴ�С��
	int RowNum = sqrt(nbrhoodD.size()-1)/2;//���ڵ�row��
	int ColNum = sqrt(nbrhoodD.size()-1)/2;//���ڵ�col��

	//double cell_res = _pDEMLayer->_pMetaData->cellSize;//	DEM��������
	double cell_res = 30;

	Pattern pattern;//pattern��ʼ��
	pattern.num_negatives = 0;
	pattern.num_positives = 0;
	pattern.negatives = NULL;
	pattern.positives = NULL;

	//���ýǶ���ֵ
	double flat_distance = cell_res*3;//Ϊ�˱��ⴰ��������������е㽫��3��֮�⿪ʼ����
	double flat_threshold = 0.01;//�������һ����
	double flat_threshold_height = cell_res*0.001;//�������һ����

	

	vector<GeoPoint> Points;//����㼯
	vector<int> PointToKeep;//���챣���㼯

	for(int i=0;i<8;i++)
	{
		pattern.pattern[i]=0;

		//��ʼ��Points
		GeoPoint curPoint;
		curPoint.i = i;

		//���ԭ��������
		Points.clear();
		PointToKeep.clear();

		for (int k = 0;k<RowNum;k++)
		{
			//cout<<"cell_res = "<<cell_res<<endl;
			curPoint.row = cur_Row + k*nextr[i];
			//cout<<"step2.1.1"<<" "<<k<<" "<<curPoint.row;
			curPoint.col = cur_Col + k*nextc[i];
			//cout<<"step2.1.2"<<" "<<k<<" "<<curPoint.col;
			curPoint.elevation = dem[curPoint.row][curPoint.col];
			//cout<<"step2.1.3"<<" "<<k<<" "<<curPoint.elevation<<" ";
			curPoint.height = curPoint.elevation - dem[cur_Row][cur_Col];
			//cout<<"step2.1.4"<<" "<<k<<" "<<curPoint.height;
			curPoint.distance = caldist(curPoint.row,curPoint.col,cur_Row,cur_Col,cell_res);
			//cout<<"step2.1 "<<curPoint.distance<<endl;
			//cout<<"curPoint.row "<<curPoint.row<<" curPoint.col "<<curPoint.col<<" cur_Row "<<cur_Row<<" cur_Col "<<cur_Col<<" cell_res "<<cell_res<<endl;
			//curPoint.distance = cell_res * sqrt((double)((curPoint.col-cur_Col)*(curPoint.col-cur_Col) +(curPoint.row-cur_Row)*(curPoint.row-cur_Row)));
			//cout<<"step2.1"<<" "<<k<<" "<<curPoint.distance<<endl;
			curPoint.angel = atan2(curPoint.height,curPoint.distance);
			Points.push_back(curPoint);
		}

		//��ʼ�����������
		double Tolerance = 0;
		int firstpoint = 0;
		int lastpoint  = Points.size()-1;

		//PointToKeep.push_back(firstpoint);
		//PointToKeep.push_back(lastpoint);

		int feature_point = DPCalcuOnetime(Points,firstpoint,lastpoint,Tolerance,PointToKeep);

		
		//cout<<feature_point<<"  "<<Points[feature_point].angel<<"  "<<flat_threshold<<endl;
		if (Points[feature_point].angel>flat_threshold)
		{
			pattern.positives += i;
			pattern.pattern[i] = 1;
			pattern.num_positives++;
			pattern.feature[i] = Points[feature_point];
		}
		if (Points[feature_point].angel<-flat_threshold)
		{
			pattern.negatives += i;
			pattern.pattern[i] = -1;
			pattern.num_negatives++;
			pattern.feature[i] = Points[feature_point];
		}
	}
	/*
	if (coord.iRow()>1)
	{
		cur_Col = coord.iRow()+100000;
	}
	*/
	//slope[cur_Row][cur_Col] = determine_form_num(pattern.num_negatives,pattern.num_positives);
	//vector<GeoPoint>().swap(Points);
	//vector<int>().swap(PointToKeep);
	slope[cur_Row][cur_Col] = determine_form_num(pattern.num_negatives,pattern.num_positives);
	return true;
}
