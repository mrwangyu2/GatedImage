# 医学影像中呼吸门控的动态显示

编写： 王宇
2018-5-25

## 呼吸门控概念
由于受到呼吸运动的影响，位于胸腹部的病灶具有较大幅度的移动。在传统放疗方式中考虑呼吸运动则需要增加较大的计划靶区（PTV），导致受照正常组织显着增加。呼吸门控技术 （respiratory gating technology）是指在放射治疗过程中监测患者呼吸运动，在特定呼吸时相触发射线照射，降低正常组织受照剂量。

## 运行效果
* 冠状面


* 矢状面



## 功能
* 解析时间、空间（x，y，z）维度数据
* 创建VTK管线，并实现动态播放
* 配置PET使用的颜色映射表（Lookup Table)
* 图片旋转180度

## Dicom相关标签
* 0054,1000 Series Type: 类型：GATED\IMAGE
* 0018,1060 Trigger Time：触发器开始获取数据的时间，时间间隔FrameTime 0
* 0018,1062 Nominal Interval: 持续时间 768
* 0018,1063 Frame Time： frame间隔时间 128
* 0018,1081 Low R-R
* 0018,1082 High R-R
* 0054,0061 Number of R-R Intervals 1
* 0054,0071 Number of Time Slots    6
* 0054,0081 Number of Slices        117
* 0054,0101 Number of Time Slices   1

## 四维数据解析原理

* 所谓的四维数据，就是时间和空间，即Trigger时间和 x, y, z
* 解析原理
  * 根据Dicom标签 Number of Time Slices 和Number of Slices 将数据进行分组， 即按照时间维度，分成不同的Series组。放在系统的临时目录temp下，然后通过vtkDICOMImageReader 读取数据文件，形成Number of Time Slots个DicomImageReader对象
  * 修改VTK中vtkIOImage项目，使其支持Dicom中关于Gated Image的标签

## 动画原理
* 解析四维数据后，系统将产生Number of Time Slots个数的Actor
* 采用了VTK InteractorStyle 的OnTime 机制，显示不同Time Slots的Actor
* 显示的方式是：从Time Slots 1 - 6 后，6 - 1的方式，这一个过程为一个呼吸


## 部署环境
* Visual Studio 2013
* 目录说明
  * data 呼吸门控的数据
  * include VTK等头文件
  * lib 编译时需要的库文件
  * dll 运行时库文件，将其复制到 proj.win32/GatedImage/Debug下
  * proj.win32 VS2013使用的工程文件
  * src 代码文件

## 程序下载地址
https://github.com/mrwangyu2/

