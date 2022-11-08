#pragma once
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cctype>
#include <sstream>
#include <algorithm>
class Preprocess
{
public:
	bool readData(std::string file_name, std::vector<std::pair<double, double>>& zeiss) {
		double x = 0;
		double y = 0;
		bool flag = false;
		std::ifstream ifs(file_name);
		if (!ifs) {
			std::cout << "File input error!" << std::endl;
			return flag;
		}
		std::string line;
		getline(ifs, line);
		while (getline(ifs, line) && line != "$END" && !line.empty()) {
			std::pair<double, double> point = transfer(line);
			zeiss.push_back(point);
		}
		Info(zeiss);
		flag = true;
		ifs.close();
		return flag;
	}

	bool Car2Pol(const std::vector < std::pair<double, double>>& car, std::vector < std::pair<double, double>>& pol) {
		for (auto point : car) {
			std::pair<double, double> postprocess(CalModule(point.first, point.second), CalAngle(point.first, point.second));
			pol.push_back(postprocess);
		}
		if (car.size() == pol.size())
			return true;
		return false;
	}

	void writeDataDev(std::string file_name, const std::vector < std::pair<double, double>>& data) {
		std::ofstream ofs(file_name);
		if (!ofs) {
			std::cout << "file open error!" << std::endl;
			return;
		}
		ofs << "$ELE,(NAM,=,PRO_M_B,TYP=APT,,FLD,=,($X,$Y,$Z))" << std::endl;
		for (const auto& point : data) {
			ofs << point.first << "," << point.second << ",0," << std::endl;
		}
		ofs << "$END" << std::endl;
	}

	void writeDataMod(std::string file_name, const std::vector < std::pair<double, double>>& data) {
		std::ofstream ofs(file_name);
		if (!ofs) {
			std::cout << "file open error!" << std::endl;
			return;
		}
		ofs << std::endl;
		for (const auto& point : data) {
			ofs << point.first << " " << point.second << " 0 " << std::endl;
		}
	}

private:

	double CalAngle(double x, double y) {
		double angle = std::atan2(y, x);
		return (angle >= 0) ? angle : 3.1415926 * 2 + angle;
	}

	double CalModule(double x, double y) {
		return sqrt(x * x + y * y);
	}

	void Info(const std::vector<std::pair<double, double>>& zeiss) {
		std::cout << "There are " << zeiss.size() << " points in the file." << std::endl;
	}

	std::pair<double, double> transfer(const std::string& str) {
		std::string num, xs, ys;
		for(auto c : str) {
			if (isdigit(c) || c == '-')
				num.push_back(c);
			if (c == ',' || isspace(c))
				num.push_back('*');
			if (c == '.')
				num.push_back('.');
		}
		auto itr = std::find(num.cbegin(), num.cend(), '*');
		for (auto i = num.cbegin(); i != itr; ++i)
			xs.push_back(*i);
		auto itr1 = std::find(++itr, num.cend(), '*');
		for (auto i = itr; i != itr1; ++i)
			ys.push_back(*i);
		double x = str2num(xs);
		double y = str2num(ys);
		return std::pair<double, double>(x, y);
	}

	double str2num(std::string str) {
		double result = 0;
		double res1 = 0;
		double res2 = 0;
		double sign = 1;
		auto b = str.cbegin();
		if (str[0] == '-') {
			b = ++(str.cbegin());
			sign = -1;
		}
		auto itr = std::find(str.cbegin(), str.cend(), '.');
		for(auto itr1 = b; itr1 != itr; ++itr1) {
			double addend = double(*itr1) - 48.0;
			res1 = res1 * 10 + addend;
		}
		if (itr != str.cend()) {
			for (auto itr2 = --str.cend(); itr2 != itr; --itr2) {
				double addend = double(*itr2) - 48.0;
				res2 = 0.1 * res2 + addend;
			}
		}
		result = sign * (res1 + 0.1 * res2);
		return result;
	}
};

