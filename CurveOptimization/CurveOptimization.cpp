﻿// CurveOptimization.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "Cycloid_Cal.h"
#include "CycloidGear.h"

using namespace std;

// 计算两条曲线的距离，并输出到文件"E:\\distance.txt"
void output_dis(vector<double> vec, std::vector<std::pair<std::pair<double, double>, double>>* flank) {
    auto itr1 = flank->cbegin();
    ofstream ofs("E:\\distance.txt");
    if (!ofs) {
        cout << "Distance Write Error!" << endl;
        return;
    }
    else {
        for (auto itr = vec.cbegin(); itr != vec.cend(); ++itr, ++itr1) {
            ofs << *itr << "," << itr1->first.first << "," << itr1->first.second << endl;
        }
    }
    return;
}

// 优化函数
void opt(const int &Z, double &R, double &r, double &e, double &obj, const double &lower, const double &upper, const string &file_name) {
    CycloidGear rv10c(Z, R / 2.0, r / 2, e);
    rv10c.setFlankRange(lower / 2, upper / 2);
    rv10c.pointInput(&file_name);

    // 摆线参数优化：
    rv10c.paramsopt();

    // 计算输入曲线与理论摆线偏差：（类似Fa）
    //output_dis(rv10c.distance(), rv10c.flank);

    // 拼接曲线（齿顶 + 优化曲线 + 齿根）：
    //rv10c.output_curve();
}

// 将数据data输出到E盘中名称为name的txt文件中
void write_data(const string& name, const vector<double>& data) {
    string path = "E:\\" + name + ".txt";
    std::ofstream ofs(path);
    if (!ofs) {
        cout << "Write Error: " << name << ".txt!" << endl;
        return;
    }
    for (auto itr = data.cbegin(); itr != data.cend(); ++itr) {
        ofs << *itr;
        ofs << endl;
    }
    return;
}

int main()
{

    const int Z = 39;
    vector<string> name{ "D121", "D122",  "D131", "D132", "D151", "D152", "D162"};                       // 点云文件名称，每个零件由一对引号标示，可放入多个零件名称。
    string dir = "E:\\WorkSpace\\02 正向研发\\摆线齿形项目\\40EPD\\";  // 点云文件绝对路径
    vector<double> R_avg = { 128.0 };                 // 针齿分布圆直径初始值
    vector<double> r_avg = { 6.0 };                   // 针齿直径初始值
    vector<double> e_avg = { 1.3 };                  // 偏心距初始值
    vector<double> obj_avg = { 0.0 };                     // 目标函数初始值
    const double lower = 120.0;                         // 评价范围下界
    const double upper = 124.0;                         // 评价范围上界
    int num = 500;                                      // 优化次数
    int n = name.size();
    for (int i = 0; i < num; i++) {
        double R_sum = 0;
        double r_sum = 0;
        double e_sum = 0;
        double obj_sum = 0;
        // 计算优化结果平均值
        for(auto str = name.begin(); str != name.end(); str++) {
            cout << "\n\niteration num: " << i << endl;
            string file_name = dir + *str + "_dev.txt";
            cout << "\n\n" << file_name << ": \n";
            double R = R_avg[i];
            double r = r_avg[i];
            double e = e_avg[i];
            double obj = obj_avg[i];
            opt(Z, R, r, e, obj, lower, upper, file_name);
            R_sum += R;
            r_sum += r;
            e_sum += e;
            obj_sum += obj;
        }
        R_avg.push_back(R_sum / n);
        r_avg.push_back(r_sum / n);
        e_avg.push_back(e_sum / n);
        obj_avg.push_back(obj_sum / n); 
        // 输出每步优化结果到文件
        write_data("R", R_avg);
        write_data("r1", r_avg);
        write_data("e", e_avg);
        write_data("obj", obj_avg);

    }

    
    
    return 0;

}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
