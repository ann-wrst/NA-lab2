#include <iostream>
#include <iomanip>
using namespace std;
const int N = 4;
const int M = 4;

struct point {
	int i;
	int j;
	double max;
};

double MatrixNorm(double** ar, int m, int n) {
	double sum;
	double max = 0;
	for (int i = 0; i < n; i++) {
		sum = 0;
		for (int j = 0; j < m; j++) {
			sum += abs(ar[i][j]);
		}
		if (sum > max)
			max = sum;
	}
	return max;
}

void RemoveFreeVariables(double ar[][N + 1], double** arFinal) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; ++j) {
			arFinal[i][j] = ar[i][j];
		}
	}
}

void MaxInColumn(int m, double ar[][N + 1], point* temp) {
	temp->max = abs(ar[0][m]);
	temp->i = 0;
	temp->j = m;
	for (int j = 0; j < N - 1; j++) {
		if (abs(ar[j + 1][m]) > temp->max) {
			temp->max = abs(ar[j + 1][m]);
			temp->i = j + 1;
			temp->j = m;
		}
	}
}

void Output(double** ar) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << ar[i][j] << ' ';
		}
		cout << endl;
	}
}
void GaussMethod(double ar[][N + 1]) {
	cout << "Gauss method" << endl;
	double** arD;
	arD = new double* [N];
	for (int i = 0; i < N; i++) {
		arD[i] = new double[N];
	}
	RemoveFreeVariables(ar, arD);
	double** e;
	e = new double* [N];
	for (int i = 0; i < N; i++) {
		e[i] = new double[N];

	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			e[i][j] = 0.0;
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			e[i][i] = 1;
		}
	}
	point Point[N];
	for (int s = 0; s < N - 1; s++) {
		MaxInColumn(s, ar, &Point[s]);
		if (Point[s].i != s) {
			for (int k = 0; k < N + 1; k++) {
				swap(ar[Point[s].i][k], ar[s][k]);
			}
			for (int k = 0; k < N; k++) {
				swap(arD[Point[s].i][k], arD[s][k]);
				swap(e[Point[s].i][k], e[s][k]);
			}
			Point[s].i = s;
		}
		for (int k = s + 1; k < N; k++) {
			if (ar[k][Point[s].j] != 0) {
				double m = -(double)ar[k][Point[s].j] / ar[Point[s].i][Point[s].j];
				for (int i = 0; i < N + 1; i++) {
					ar[k][i] += ar[Point[s].i][i] * m;
				}
				for (int i = 0; i < N; i++) {
					e[k][i] += e[Point[s].i][i] * m;
					arD[k][i] += arD[Point[s].i][i] * m;
				}
			}
		}
	}
	double det = 1;
	double dg;
	for (int i = 0; i < N; i++) {
		dg = arD[i][i];
		det *= dg;
	}
	cout << fixed << setprecision(4);
	cout << "The determinant of the matrix = " << det << endl << endl;

	for (int i = 0; i < N; i++) {
		double diag = arD[i][i];
		for (int j = 0; j < N; j++) {
			if (diag != 0) {
				arD[i][j] /= diag;
				e[i][j] /= diag;
			}
		}
	}

	for (int j = N - 1; j >= 1; j--) {
		for (int i = j - 1; i >= 0; i--) {
			double m = -arD[i][j];
			for (int k = 0; k < N; k++) {
				arD[i][k] += arD[j][k] * m;
				e[i][k] += e[j][k] * m;
			}
		}
	}
	cout << "The inverse matrix is\n";
	cout << endl;
	Output(e);
	cout << endl;
	double conditionNumber = MatrixNorm(arD, N, N) * MatrixNorm(e, N, N);
	cout << "The condition number is " << conditionNumber << endl << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N + 1; j++) {
			cout << ar[i][j] << ' ';
		}
		cout << endl;
	}
	double x[N];
	x[3] = ar[N - 1][N] / ar[N - 1][N - 1];
	x[2] = (ar[N - 2][N] - x[3] * ar[N - 2][N - 1]) / ar[N - 2][N - 2];
	x[1] = (ar[N - 3][N] - x[3] * ar[N - 3][N - 1] - x[2] * ar[N - 3][N - 2]) / ar[N - 3][N - 3];
	x[0] = (ar[N - 4][N] - x[3] * ar[N - 4][N - 1] - x[2] * ar[N - 4][N - 2] - x[1] * ar[N - 4][N - 3]) / (double)(ar[N - 4][N - 4]);
	cout << endl;
	for (int i = 0; i < N; i++) {
		cout << "x" << i + 1 << " = " << x[i] << " ";
	}
	cout << endl;
}

void JacobiMethod(double ar[][M + 1], double x[], double xR[], double q, double eps = 1e-3) {
	cout << "\nJacobi method\n";
	for (int i = 0; i < M; i++) {
		if (ar[i][i] == 0) {
			cout << "Diagonal element = 0" << endl;
			return;
		}
	}
	for (int i = 0; i < M; i++) {
		double sum = 0;
		double currentDiagEl = ar[i][i];
		for (int j = 0; j < M; j++) {
			sum += abs(ar[i][j]);
			if (currentDiagEl < sum - currentDiagEl) {
				cout << "The matrix is not diagonally dominant" << endl;
				return;
			}
		}
	}
	double div;
	for (int i = 0; i < M; i++) {
		div = ar[i][i];
		for (int j = 0; j < M + 1; j++) {
			ar[i][j] /= div;
		}
	}
	int n;
	int choice;
	double x_res[4];
	cout << "Do you want to enter the number of iterations?\n1.Yes\n2.No" << endl;
	cin >> choice;
	if (choice == 1) {
		cout << "Enter the number of iterations" << endl;
		cin >> n;
		int k = 0;
		cout << "n\t";
		for (int i = 0; i < M; i++) {
			cout << "x" << i + 1 << "\t";
		}
		cout << endl;
		while (k < n) {
			cout << endl << k + 1 << "\t";
			for (int i = 0; i < M; i++) {
				std::pair<double, double> sum;
				sum.first = 0;
				sum.second = 0;

				for (int j = 0; j < i; j++)
					sum.first += ar[i][j] * x[j];

				for (int j = i + 1; j < M; j++)
					sum.second += ar[i][j] * x[j];

				x_res[i] = -sum.first - sum.second + ar[i][4];
				cout << fixed;
				cout << setprecision(4) << x_res[i] << "\t";
			}
			cout << endl;
			x[0] = x_res[0];
			x[1] = x_res[1];
			x[2] = x_res[2];
			x[3] = x_res[3];
			k++;
		}
	}
	else
		if (choice == 2) {
			cout << "n\t";
			for (int i = 0; i < M; i++) {
				cout << "x" << i + 1 << "\t";
			}
			cout << endl;
			int count = 0;
			for (int m = 0; m < M; m++) {
				while (abs(x_res[m] - xR[m]) >= eps) {
					count++;
					cout << endl << count << "\t";
					for (int i = 0; i < M; i++) {
						std::pair<double, double> sum;
						sum.first = 0;
						sum.second = 0;

						for (int j = 0; j < i; j++)
							sum.first += ar[i][j] * x[j];

						for (int j = i + 1; j < M; j++)
							sum.second += ar[i][j] * x[j];

						x_res[i] = -sum.first - sum.second + ar[i][4];
						cout << fixed;
						cout << setprecision(4) << x_res[i] << "\t";
					}
					cout << endl;
					x[0] = x_res[0];
					x[1] = x_res[1];
					x[2] = x_res[2];
					x[3] = x_res[3];
				}
			}

		}
		else cout << "Incorrect input" << endl;
}
int main()
{
	double eqGauss[N][N + 1] = {
	{7.0,2.0,3.0,0.0,32.0},
	{0.0,3.0,2.0,6.0,47.0},
	{2.0,5.0,1.0,0.0,23.0},
	{0.0,1.0,4.0,2.0,29.0}
	};
	double eqJacobi[M][M + 1] = {
	{4.0,0.0,1.0,0.0,12.0},
	{0.0,3.0,0.0,2.0,19.0},
	{1.0,0.0,5.0,1.0,27.0},
	{0.0,2.0,1.0,4.0,30.0}
	};
	double x0[4] = { 0,0,0,0 };
	double xR[4] = { 2,3,4,5 };
	GaussMethod(eqGauss);
	JacobiMethod(eqJacobi, x0, xR, 0.75);
	system("pause");
}