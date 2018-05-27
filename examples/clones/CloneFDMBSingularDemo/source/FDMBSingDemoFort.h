//
// 文件名: FDMBSingDemoFort.h
// 软件包: JASMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1.2 $
// 修改  : $Date: 2007/01/06 09:07:25 $
// 描述  : F77外部子程序接口.
//

#include <math.h>
#include <signal.h>

extern "C" {
 
  void summing_(
  const int& , const int& , 
  const int& , const int& , 
  const int& , const int& ,  
  const double*,
  double*);

}
