# DPGeomorphons
Geomorphons算法使用并行算子PaRGO实现，但是其中寻找特征点的方法并非使用“地形开放度”算法，而是使用“Douglas-Peucker”算法。   

----使用DP算法的原始Geomorphons算法----   
----王彦文，2015年----   

##重要文件    
###1. pSlope.cpp是主函数，含有算法入口。（该程序是根据占立军师兄计算坡度算法改写，因此主程序名称为pSlope.cpp）。    
###2. slopeOperator.cpp与slopeOperater.h是判断栅格属于何种地形元素的函数，geomorphons.h是上一过程需要用到的函数集。
###3. neigh.nbr是窗口大小文件，例如半径为60的窗口，则值为(2*60+1)^2=14641
