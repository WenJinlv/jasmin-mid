# JASMIN 精简版

JASMIN-Lite 是 JASMIN 框架的单层单块网格精简版本。它保留了 JASMIN 框架的基础架 构，移除了多块、自适应、联邦、克隆、扫描和解法器功能。

## 功能特性

JASMIN-Lite 基于 MPI+OpenMP 混合编程模型，支持大规模科学工程计算应用。 JASMIN-Lite 实现了如下的特性：

1. 支持结构网格构件化编程

  JASMIN-Lite 将结构网格计算中共性的并行技术抽象成 `积分构件`，用户只需要实 现构件要求的串行策略函数，组合使用构件实现并行结构网格计算。JASMIM-Lite 包 含的构件有：

  - NumericalIntegratorComponent：数值构件
  - ReductionIntegratorComponent：规约构件
  - InitializeIntegratorComponent：初值构件
  - DtIntegratorComponent：步长构件
  - MemoryIntegratorComponent：内存构件
  - CopyIntegratorComponent：复制构件
  - OuterdataOperationIntegratorComponent: 外表面数据操作构件
  - ParticleCommComponent：粒子量通信构件
  - RemappingIntegratorComponent：重分重映构件
  - LocalNumericalIntegratorComponent：局部数值构件
  - LocalReductionIntegratorComponent：局部规约构件
  - StencilNumericalIntegratorComponent：Stencil 模板数值构件

2. 支持多种结构网格应用开发

  JASMIN-Lite 支持如下应用：

  - 单块均匀矩形网格应用
  - 单块变形网格应用和网格重分重映
  - PIC方法和粒子模拟应用

3. 支持动态负载平衡

## 编译安装

从源代码编译 JASMIN-Lite 需要系统已经安装好如下软件包：

1. MPI & OpenMP：建议安装 mpich2 或 openmpi；建议 gcc 版本 4.4 以上；编译器 必须支持OpenMP 3.0 标准。
2. HDF5 1.6 版本
3. FFTW 3.0 以上版本
4. BOOST 1.44 以上版本

在 Linux 环境中，通过如下的命令编译安装 JASMIN-Lite：

```bash
mkdir build && cd build
cp $SOURCE/tools/config-linux-optg.sh .
# 修改 config-linux-optg.sh 中的文件路径
./config-linux-optg.sh
make && make install
```

## 使用方法

使用方法参见安装目录中 `examples/*`。

## 功能详解

1. Toolbox

  提供基本的辅助功能支持，对用户提供如下的能力：

  1. 科学计算数据输出: hdf5, jdio
  2. 断点续算功能: restart database
  3. (Optional) 和作业系统配合，容错计算功能：应用程序长时间无干预持续计算 计算。
  4. 基本的性能监测和热点分析功能
  5. 通过输入文件描述计算模型的能力
