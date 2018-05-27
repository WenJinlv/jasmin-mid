//
// File:        DynamicFort.h
// Package:     JASMIN application
// Copyright:   (c) 2009-1014 The Regents of IAPCM
// Revision:    $Revision: 1.8 $
// Modified:    $Date: 2007/08/20 01:39:33 $
// Description: 声明求解Dynamic模型问题的F77子程序.
//

extern "C"
  {
    // 非协调边界垂直于X方向的问题, 初始化网格片数据, 计算区域包含两个block.
    void initializefield2blockv_(
      const int&, const int&,      // ifirst(0),ilast(0)
      const int&, const int&,      // ifirst(1),ilast(1)
      const int&, const int&,      // ghost_cells(0),ghost_cells(1)
      const int&,                  // block_id
      const double*,               // d_velocity_in
      double*,                     // velocity
      double*);                    // node_coord

    // 非协调边界垂直于Y方向的问题, 初始化网格片数据, 计算区域包含两个block.
    void initializefield2blockh_(
      const int&, const int&,      // ifirst(0),ilast(0)
      const int&, const int&,      // ifirst(1),ilast(1)
      const int&, const int&,      // ghost_cells(0),ghost_cells(1)
      const int&,                  // block_id
      const double*,               // d_velocity_in
      double*,                     // velocity
      double*);                    // node_coord

    // 非协调边界垂直于Y方向的问题, 初始化网格片数据, 计算区域包含三个block.
    void initializefield3blockh_(
      const int&, const int&,      // ifirst(0),ilast(0)
      const int&, const int&,      // ifirst(1),ilast(1)
      const int&, const int&,      // ghost_cells(0),ghost_cells(1)
      const int&,                  // block_id
      const double*,               // d_velocity_in
      double*,                     // velocity
      double*);                    // node_coord

    // 计算新时刻的结点位置坐标.
    void advancecoordinate_(
      const int&, const int&,      // ifirst(0),ilast(0)
      const int&, const int&,      // ifirst(1),ilast(1)
      const int&, const int&,      // ghost_cells(0),ghost_cells(1)
      const double&,               // dt
      const double*,               // velocity
      double*);                    // node_coord
  }
