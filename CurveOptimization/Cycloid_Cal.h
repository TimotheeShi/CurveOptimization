#ifndef CYCLOID_CAL_H
#define CYCLOID_CAL_H

#include <utility>
#include <vector>
#include <math.h>

/**
 * 此头文件用于计算理论摆线及其导数
 */


/**
 * 该函数为理论摆线的参数方程
 * @params:
 *		t: 参数方程的相位角；
 *		R: 针齿分布圆直径；
 *		Z: 摆线齿轮齿数；
 *		r: 针齿直径；
 *		e: 偏心距。
 * @returns:
 *		返回摆线上点的极坐标（row，theta）。
 */
std::pair<double, double> cycloid_pol(const double& t, const double& R, const double& Z, const double& r, const double& e) {
	double x = 0.0, y = 0.0;
	double Z2 = Z + 1;
	double K = e * Z2 / R;
	double S = 1 + K * K - 2 * K * cos(Z * t);
	double sgamma = (-K * cos(Z2 * t) + cos(t)) / sqrt(S);
	double cgamma = (K * sin(Z2 * t) - sin(t)) / sqrt(S);
	double x0 = R * cos(t) - e * cos(Z2 * t);
	double y0 = R * sin(t) - e * sin(Z2 * t);
	x = x0 - r * sgamma;
	y = y0 + r * cgamma;
	double row = sqrt(x * x + y * y);
	return std::pair<double, double> (row, t);
}

/**
 * 该函数为理论摆线的参数方程
 * @params:
 *		t: 参数方程的相位角；
 *		R: 针齿分布圆直径；
 *		Z: 摆线齿轮齿数；
 *		r: 针齿直径；
 *		e: 偏心距。
 * @returns:
 *		返回摆线上点的笛卡尔坐标: (x, y)。
 */
std::pair<double, double> cycloid_car(const double& t, const double& R, const double& Z, const double& r, const double& e) {
	double x = 0.0, y = 0.0;
	double Z2 = Z + 1;
	double K = e * Z2 / R;
	double S = 1 + K * K - 2 * K * cos(Z * t);
	double sgamma = (-K * cos(Z2 * t) + cos(t)) / sqrt(S);
	double cgamma = (K * sin(Z2 * t) - sin(t)) / sqrt(S);
	double x0 = R * cos(t) - e * cos(Z2 * t);
	double y0 = R * sin(t) - e * sin(Z2 * t);
	x = x0 - r * sgamma;
	y = y0 + r * cgamma;
	return std::pair<double, double>(x, y);
}


/**
 * 该函数为计算了摆线的模分别对R,r,e求导的导数值
 * @params:
 *		phi: 参数方程的相位角；
 *		R: 针齿分布圆直径；
 *		Z: 摆线齿轮齿数；
 *		r: 针齿直径；
 *		e: 偏心距。
 * @returns:
 *		返回一个大小为3×1的Jacobian矩阵：[dρ/dR, dρ/dr, dρ/de]'。(ρ为模)
 */
std::vector<double> cycloid_pol_derivative(const double& t, const double& R, const double& Z, const double& r, const double& e) {
	double dR = 0.0, dr = 0.0, de = 0.0;
	double x = 0.0, y = 0.0;
	double Z2 = Z + 1;
	double K = e * Z2 / R;
	double S = 1 + K * K - 2 * K * cos(Z * t);
	double sgamma = (-K * cos(Z2 * t) + cos(t)) / sqrt(S);
	double cgamma = (K * sin(Z2 * t) - sin(t)) / sqrt(S);
	double x0 = R * cos(t) - e * cos(Z2 * t);
	double y0 = R * sin(t) - e * sin(Z2 * t);
	x = x0 - r * sgamma;
	y = y0 + r * cgamma;
	double row = sqrt(x * x + y * y);

	// Calculate drow/dr
	double drow_dx = x / row;
	double drow_dy = y / row;
	double dx_dr = -sgamma;
	double dy_dr = cgamma;
	dr = drow_dx * dx_dr + drow_dy * dy_dr;

	// Calculate drow/dR
	double dx_dx0 = 1, dy_dy0 = 1;
	double dx_dsgamma = -r;
	double dy_dcgamma = r;
	double dx0_dR = cos(t);
	double dy0_dR = sin(t);
	double dK_dR = -e * Z2 / (R * R);
	double dS_dR = 2 * K * dK_dR - 2 * dK_dR * cos(Z * t);
	double dsgamma_dR = (-dK_dR * cos(Z2 * t) * sqrt(S) - (-K * cos(Z2 * t) + cos(t)) * (0.5 / sqrt(S)) * dS_dR) / S;
	double dcgamma_dR = (dK_dR * sin(Z2 * t) * sqrt(S) - (K * sin(Z2 * t) - sin(t)) * (0.5 / sqrt(S)) * dS_dR) / S;
	dR = drow_dx * (dx_dx0 * dx0_dR + dx_dsgamma * dsgamma_dR) + drow_dy * (dy_dy0 * dy0_dR + dy_dcgamma * dcgamma_dR);

	// Calculate drow/de
	double dK_de = Z2 / R;
	double dS_de = 2 * K * dK_de - 2 * cos(Z * t) * dK_de;
	double dx0_de = -cos(Z2 * t);
	double dy0_de = -sin(Z2 * t);
	double dsgamma_de = (-dK_de * cos(Z2 * t) * sqrt(S) - (-K * cos(Z2 * t) + cos(t)) * 0.5 / sqrt(S) * dS_de) / S;
	double dcgamma_de = (dK_de * sin(Z2 * t) * sqrt(S) - (K * sin(Z2 * t) - sin(t)) * 0.5 / sqrt(S) * dS_de) / S;
	de = drow_dx * (dx_dx0 * dx0_de + dx_dsgamma * dsgamma_de) + drow_dy * (dy_dy0 * dy0_de + dy_dcgamma * dcgamma_de);
	
	return std::vector<double> {dR, dr, de};
}

/**
 * 该函数为计算了摆线的模分别对R,r,e求导的导数值
 * @params:
 *		phi: 参数方程的相位角；
 *		R: 针齿分布圆直径；
 *		Z: 摆线齿轮齿数；
 *		r: 针齿直径；
 *		e: 偏心距。
 * @returns:
 *		返回一个大小为3×2的Jacobian矩阵：[dx/dR, dy/dR; dx/dr, dy/dr; dx/de, dy/de]。
 */
std::vector<std::pair<double, double>> cycloid_car_derivative(const double& t, const double& R, const double& Z, const double& r, const double& e) {
	std::vector<std::pair<double, double>> derivative;
	
	double dR = 0.0, dr = 0.0, de = 0.0;
	double x = 0.0, y = 0.0;
	double Z2 = Z + 1;
	double K = e * Z2 / R;
	double S = 1 + K * K - 2 * K * cos(Z * t);
	double sgamma = (-K * cos(Z2 * t) + cos(t)) / sqrt(S);
	double cgamma = (K * sin(Z2 * t) - sin(t)) / sqrt(S);
	double x0 = R * cos(t) - e * cos(Z2 * t);
	double y0 = R * sin(t) - e * sin(Z2 * t);
	x = x0 - r * sgamma;
	y = y0 + r * cgamma;

	// Calculate dx/dr, dy/dr
	double dx_dr = -sgamma;
	double dy_dr = cgamma;

	// Calculate dx/dR, dy/dR
	double dx_dx0 = 1, dy_dy0 = 1;
	double dx_dsgamma = -r;
	double dy_dcgamma = r;
	double dx0_dR = cos(t);
	double dy0_dR = sin(t);
	double dK_dR = -e * Z2 / (R * R);
	double dS_dR = 2 * K * dK_dR - 2 * dK_dR * cos(Z * t);
	double dsgamma_dR = (-dK_dR * cos(Z2 * t) * sqrt(S) - (-K * cos(Z2 * t) + cos(t)) * (0.5 / sqrt(S)) * dS_dR) / S;
	double dcgamma_dR = (dK_dR * sin(Z2 * t) * sqrt(S) - (K * sin(Z2 * t) - sin(t)) * (0.5 / sqrt(S)) * dS_dR) / S;
	double dx_dR = dx_dx0 * dx0_dR + dx_dsgamma * dsgamma_dR;
	double dy_dR = dy_dy0 * dy0_dR + dy_dcgamma * dcgamma_dR;

	// Calculate dx/de, dy/de
	double dK_de = Z2 / R;
	double dS_de = 2 * K * dK_de - 2 * cos(Z * t) * dK_de;
	double dx0_de = -cos(Z2 * t);
	double dy0_de = -sin(Z2 * t);
	double dsgamma_de = (-dK_de * cos(Z2 * t) * sqrt(S) - (-K * cos(Z2 * t) + cos(t)) * 0.5 / sqrt(S) * dS_de) / S;
	double dcgamma_de = (dK_de * sin(Z2 * t) * sqrt(S) - (K * sin(Z2 * t) - sin(t)) * 0.5 / sqrt(S) * dS_de) / S;
	double dx_de = dx_dx0 * dx0_de + dx_dsgamma * dsgamma_de;
	double dy_de = dy_dy0 * dy0_de + dy_dcgamma * dcgamma_de;

	// Save derivative results
	derivative.push_back(std::pair<double, double>(dx_dR, dy_dR));
	derivative.push_back(std::pair<double, double>(dx_dr, dy_dr));
	derivative.push_back(std::pair<double, double>(dx_de, dy_de));

	return derivative;
}

#endif // !CYCLOID_CAL.H
