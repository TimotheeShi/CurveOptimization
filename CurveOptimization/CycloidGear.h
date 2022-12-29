#pragma once
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include "Preprocess.h"
#include "PointInspection.h"
#include "Cycloid_Cal.h"
#include <map>

using namespace std;

class CycloidGear
{
public:
	// ���캯���������ڳݲ���ʼ����������
	CycloidGear(double Z1, double R1, double r1, double e1) : R(R1), R0(R1), Z(Z1), r(r1), r0(r1), e(e1), e0(e1) {
		zeiss = new std::vector<std::pair<double, double>>;
		zeiss_pol = new std::vector<std::pair<double, double>>;
		zeiss_opt = new std::vector<std::pair<std::pair<double, double>, double>>;
		flank = new std::vector<std::pair<std::pair<double, double>, double>>;
		addendum = new std::vector<std::pair<std::pair<double, double>, double>>;
		root = new std::vector<std::pair<std::pair<double, double>, double>>;
	};

	// �����������ͷ���ռ�õ��ڴ���Դ
	~CycloidGear() {
		delete zeiss;
		delete zeiss_pol;
		delete flank;
		delete addendum;
		delete root;
		zeiss = nullptr;
		zeiss_pol = nullptr;
		flank = nullptr;
		addendum = nullptr;
		root = nullptr;
	}

	// �趨�ݲ෶Χ�������۷�Χ
	void setFlankRange(const double lower, const double upper) {
		this->lower = lower;
		this->upper = upper;
	}

	// ���ƶ�ȡ����
	void pointInput(const string* path = nullptr) {
		if (path == nullptr) {
			cout << "�յ�ַ��" << endl;
			return;
		}
		Preprocess pp;
		PointInspection pi;
		if (pp.readData(*path, *zeiss)) {
			if (pp.Car2Pol(*zeiss, *zeiss_pol)) {
				cout << "Read Successfully!" << endl;
				theoritical_ang(zeiss, zeiss_opt);
				for (int i = 0; i < zeiss_pol->size(); i++) {
					if ((*zeiss_pol)[i].first >= this->lower && (*zeiss_pol)[i].first <= this->upper) {
						flank->push_back((*zeiss_opt)[i]);
					}
					if ((*zeiss_pol)[i].first > this->upper) {
						addendum->push_back((*zeiss_opt)[i]);
					}
					if ((*zeiss_pol)[i].first < this->lower) {
						root->push_back((*zeiss_opt)[i]);
					}
				}
			}
		}
		else
			cout << "Fail to Read!" << endl;
	}

	// ���� sqrt(x * x + y * y)
	double mod(const double x, const double y) {
		return sqrt(x * x + y * y);
	}

	// �����Ƶĵ�˳����������Ż��㷨
	void shuffle(std::vector<std::pair<double, double>>* points) {
		std::vector<std::pair<double, double>>* new_points = new std::vector<std::pair<double, double>>;
		for (int i = 0; i < Z; i++) {
			for (int j = i; j < points->size(); j += Z) {
				new_points->push_back((*points)[j]);
			}
		}
		std::swap(new_points, points);
	}

	void paramsOutput() {
		cout << "���߳��ֳ���Ϊ��" << Z;
		cout << "\n���۰��߲���Ϊ��\n\tR = " << 2*R0 << "\n\tr = " << 2*r0 << "\n\te = " << e0;
		cout << "\nʵ�ʰ��߲���Ϊ��\n\tR = " << 2*R << "\n\tr = " << 2*r << "\n\te = " << e;
		cout << endl;
	}

	// ��������ϵĵ�����ۼ���
	bool theoritical_ang(std::vector<std::pair<double, double>>* zeiss, std::vector<std::pair<std::pair<double, double>, double>>* zeiss_opt) {
		for (const pair<double, double>& point : *zeiss) {
			double ang0 = atan2(point.second, point.first);
			int acc = 3; // According to experiments, acc must smaller than 3!
			double range = 3.1415926 / Z / acc; // Accelerate the finding process. Higher the acc is, faster the process is.
			double ang = ang0 - range ;
			double min_dis = distance(point.first, point.second, 0, 0);
			double min_ang = ang;
			for (; ang < ang0 + range ; ang += 0.000001) { // 0.000001
				auto the = cycloid_car(ang, R, Z, r, e);
				double curr_dis = distance(point.first, point.second, the.first, the.second);
				if (min_dis > curr_dis) {
					min_dis = curr_dis;
					min_ang = ang;
				}
			}
			zeiss_opt->push_back(std::pair<std::pair<double, double>, double>(point, min_ang));
		}
		if (zeiss->size() == zeiss_opt->size())
			return true;
		return false;
	}

	// ��������(x1, y1) & (x2, y2)����
	double distance(const double x1, const double y1, const double x2, const double y2) {
		double dx = x1 - x2;
		double dy = y1 - y2;
		return sqrt(dx * dx + dy * dy);
	}

	// �Ż�����
	void paramsopt() {
		double obj = 10000;
		double old_obj = 0;
		double diff_obj = 1000;
		int j = 0;
		int k = 0;
		while (diff_obj >= 1e-10) {
			double old_dR = 0;
			double old_de = 0;
			double old_dr = 0;
			vector<double> res_x_squre, res_y_squre;
			vector<double> jacobian_R;
			vector<double> jacobian_r;
			vector<double> jacobian_e;

			// choose the data set
			auto points = flank;

			// define the residual and the object function
			for (auto pt = points->cbegin(); pt != points->cend(); ++pt) {
				pair<double, double> theo = cycloid_car(pt->second, this->R, this->Z, this->r, this->e);
				double res_x = pt->first.first - theo.first;
				double res_y = pt->first.second - theo.second;
				res_x_squre.push_back(res_x * res_x);
				res_y_squre.push_back(res_y * res_y);
				auto derivative = cycloid_car_derivative(pt->second, this->R, this->Z, this->r, this->e);
				jacobian_R.push_back(-res_x * derivative[0].first - res_y * derivative[0].second);
				jacobian_r.push_back(-res_x * derivative[1].first - res_y * derivative[1].second);
				jacobian_e.push_back(-res_x * derivative[2].first - res_y * derivative[2].second);
			}
 			old_obj = obj;
			obj = 0.5 * vec_sum(res_x_squre) + 0.5 * vec_sum(res_y_squre);

			// update optimization parameters
			double dR = vec_sum(jacobian_R);
			double dr = vec_sum(jacobian_r);
			double de = vec_sum(jacobian_e);
			double lr = -0.00001;
			R += lr * dR;
			r += lr * dr;
			e += lr * de;

			diff_obj = old_obj - obj;
			cout << "[" << ++j << "] Cost function: " << obj << ", diff_cost: " << diff_obj << ", dR=" << 2*R << ", dr=" << 2*r << ", e=" << e << endl;
			if (dR * old_dR < 0 || dr * old_dr < 0 || de * old_de < 0)
				break;
			old_dR = dR;
			old_dr = dr;
			old_de = de;
		}
		obj_func = obj;
	}

	// Ѱ��vector������Ԫ��
	double find_max(vector<double> vec) {
		int res = 0;
		for (auto itr = vec.cbegin(); itr != vec.cend(); ++itr) {
			if (*itr > res)
				res = *itr;
		}
		return res;
	}

	// ������Ƴ�����Χ�ڵ��������۵�֮��ľ���
	vector<double> distance() {
		vector<double> dis;
		double max = 0;
		double min = 100;
		for (auto itr = flank->cbegin(); itr != flank->cend(); itr++) {
			double d = cal_dis(itr->first, itr->second);
			dis.push_back(d);
			if (d > max)
				max = d;
			else if (d < min)
				min = d;
		}
		return dis;
	}

	// ���������ĳ���������۵����
	double cal_dis(pair<double, double> point, double ang) {
		pair<double, double> point_t = cycloid_car(ang, R, Z, r, e);
		int flag = 1;
		double d_x = point.first - point_t.first;
		double d_y = point.second - point_t.second;
		if (sqrt(point.first * point.first + point.second * point.second) - sqrt(point_t.first * point_t.first + point_t.second * point_t.second) < 0)
			flag = -1;
		return flag * sqrt(d_x * d_x + d_y * d_y);
	}

	// ������Ƶ���������������
	pair<double, double> cal_modification() {
		std::pair<double, double> addendum = (*zeiss_pol)[0];
		std::pair<double, double> root = (*zeiss_pol)[0];
		for (auto itr = zeiss_pol->cbegin(); itr != zeiss_pol->cend(); ++itr) {
			if (itr->first > addendum.first)
				addendum = *itr;
			if (itr->first < root.first)
				root = *itr;
		}
		double add_mod = cycloid_pol(addendum.second, this->R, this->Z, this->r, this->e).first - addendum.first;
		double root_mod = cycloid_pol(root.second, this->R, this->Z, this->r, this->e).first - root.first;
		return pair<double, double>(add_mod, root_mod);
	}

	// ��vector������Ԫ�صĴ�����
	double vec_sum(const vector<double>& vec) {
		double sum = 0.0;
		for (auto itr = vec.cbegin(); itr != vec.cend(); ++itr)
			sum += *itr;
		return sum;
	}

	// ���ƴ�����ߣ��ݶ� + �ݲ� + �ݸ���
	void output_curve() {
		std::map<double, double> output;
		for (auto itr = addendum->cbegin(); itr != addendum->cend(); ++itr) {
			if (itr->second <= 3.14159 / 39 + 3.14159 / 180 && itr->second >= - 3.14159 / 39 - 3.14159 / 180)
				output[itr->first.second] = itr->first.first;
		}
		for (auto itr = root->cbegin(); itr != root->cend(); ++itr) {
			if (itr->second <= 3.14159 / 39 + 3.14159 / 180 && itr->second >= -3.14159 / 39 - 3.14159 / 180)
				output[itr->first.second] = itr->first.first;
		}
		for (auto itr = flank->cbegin(); itr != flank->cend(); ++itr) {
			if (itr->second <= 3.14159 / 39 + 3.14159 / 180 && itr->second >= -3.14159 / 39 - 3.14159 / 180)
				output[cycloid_car(itr->second, R, Z, r, e).second] = cycloid_car(itr->second, R, Z, r, e).first;
		}

		std::ofstream ofs("E:\\output.txt");
		ofs << "$ELE (NAM = CIR_11(1), TYP=APT, FLD = ($X, $Y))\n";
		for (auto itr = output.cbegin(); itr != output.cend(); ++itr) {
			ofs << itr->second << ", " << itr->first << ", 0" << std::endl;
		}
		ofs << "$END";
		ofs.close();

	}

	// variables with the postfix "0" are the theoretical values (which are also initial values of the optimition).
	double R = 0, R0 = 0;
	double Z = 0;
	double r = 0, r0 = 0;
	double e = 0, e0 = 0;
	double upper = R;
	double lower = 0;
	double obj_func = 0;
	std::vector<std::pair<double, double>>* zeiss;
	std::vector<std::pair<double, double>>* zeiss_pol;
	std::vector<std::pair<std::pair<double, double>, double>>* zeiss_opt;
	std::vector<std::pair<std::pair<double, double>, double>>* flank;
	std::vector<std::pair<std::pair<double, double>, double>>* addendum;
	std::vector<std::pair<std::pair<double, double>, double>>* root;
	std::vector<std::pair<std::pair<double, double>, double>>* output_flank;

	// ��ȡd�ķ���
	double signs(double d) {
		double res = 0;
		if (d >= 0)
			res = 1;
		if (d < 0)
			res = -1;
		return res;
	}

};

