#ifndef SLOPEOPERATOR_H
#define SLOPEOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class SlopeOperator : public RasterOperator<double> 
{
  public:
    SlopeOperator()
      :RasterOperator<double>(),
       _pDEMLayer(0), _pSlopeLayer(0), num(0){}
   
    ~SlopeOperator() {}

  
    void demLayer(RasterLayer<double> &layerD);
	void slopeLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();

    virtual bool Operator(const CellCoord &coord);

  protected:
	int cellSize;
	int noData;
	int num;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pSlopeLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif