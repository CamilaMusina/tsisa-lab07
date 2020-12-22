#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <random>
#include <stdexcept>

struct Value
{
    double h;
    double dist;
    std::vector<double> al;
    double w;
    double d;
};

double n(const double& P, const double& E, const double& Xmax, const double& Xmin)
{
    return (log(1 - P) / log(1 - E / (Xmax - Xmin)));
}

std::vector<double> alpha(const int& r)
{
    std::vector<double> al(r);
    std::mt19937 gen(time(0));
    std::uniform_real_distribution<double> urd(0, 1);
    double A = urd(gen);
    al[(r - 1) / 2] = A;

    for (auto i = (r - 1) / 2 - 1; i > 0; i--)
    {
        double sum = 0;
        for (auto j = i; j < r - i - 1; j++)
            sum += al[j];
        std::uniform_real_distribution<double> urd1(0, 1 - sum);
        double B = urd1(gen);
        al[i] = 0.5 * B;
    }
    for (auto i = (r - 1) / 2 + 1; i < r - 1; i++)
        al[i] = al[r - i - 1];
    double sum = 0;
    for (auto i = 1; i < r - 1; i++)
        sum += al[i];
    al[0] = 0.5 * (1 - sum);
    al[r - 1] = al[0];
    return al;
}

std::vector<double> filter(const std::vector<double>& F, std::vector<double>& al, const int& r)
{
    std::vector<double> F_res(101);
    int M = (r - 1) / 2;
    for (auto i = 0; i <= 100; i++)
    {
        double sum = 0;
        for (auto j = i - M; j <= i + M; j++)
        {
            try
            {
                F.at(j);
                sum += F[j] * al[j + M - i];
            }
            catch (const std::out_of_range& s)
            {
                sum += 0;
            }
        }
        F_res[i] = sum;
    }
    return F_res;
}

double w(const std::vector<double>& F)
{
    double sum = 0;
    for (auto i = 1; i <= 100; i++)
        sum += pow((F[i] - F[i - 1]), 2);
    return sqrt(sum);
}

double d(const std::vector<double>& F, const std::vector<double>& F_noise)
{
    double sum = 0;
    for (auto i = 1; i <= 100; i++)
        sum += pow((F[i] - F[i - 1]), 2);
    return sqrt(sum/100);
}

double dist(const double& w, const double& d)
{
    return sqrt(pow(w, 2) + pow(d, 2));
}

double J(double h, double w, double d)
{
    return h * w + (1 - h) * d;
}

int main()
{
    int r3 = 3, r5 = 5;
    double Xmin = 0, Xmax = M_PI;
    int K = 100;
    double a1 = -0.25, a2 = 0.25;
    int L = 10;
    double H = 0;
    double P = 0.95;
    double E = 0.01;
    double W, D, minW, minD;
    std::vector<Value> values;
    std::vector<double> al;
    std::vector<double> almin;
    double Dist, minDist = 100000;
    double N = n(P, E, Xmax, Xmin);

    std::vector<double> X(101);
    for (auto i = 0; i < K; i++)
        X[i] = Xmin + i * (Xmax - Xmin) / K;

    std::vector<double> F(101);
    for (auto i = 0; i <= K; i++)
        F[i] = sin(X[i]) + 0.5;

    std::vector<double> F_noise(101);
    for (auto i = 0; i <= K; i++)
    {
        std::mt19937 gen(time(0));
        std::uniform_real_distribution<double> urd(a1, a2);
        double A = urd(gen);
        F_noise[i] = F[i] + A;
    }

    for (auto i = 0; i <= K; i++)
        std::cout << "F = " << F[i] << " F with noise = " << F_noise[i] << "\n";


    std::vector<double> F_r3;
    std::cout << " h    dist             alpha              w        d" << "\n";
    for (auto l = 0; l <= L; l++)
    {
        for (auto i = 0; i < K / 2; i++)
        {
            al = alpha(r3);
            F_r3 = filter(F_noise, al, r3);
            W = w(F_r3);
            D = d(F_r3, F_noise);
            Dist = dist(W, D);

            if (Dist < minDist)
            {
                minDist = Dist;
                minW = W;
                minD = D;
                almin = al;
            }
        }
        H = (double)l / L;
        values.push_back({ H, minDist, almin, minW, minD });
        std::cout << std::fixed << std::setprecision(1) << H << "  " << std::setprecision(4) << minDist << "   [ ";
        for (auto i : almin) 
            std::cout << i << " ";
        std::cout << "]  " << minW << "   " << minD << "\n";
        minDist = 100000;
    } 

    std::sort(values.begin(), values.end(), [](auto x, auto y) {return x.dist < y.dist; });
    std::cout << "\n" << " h*     J        w        d" << "\n";
    std::cout << std::setprecision(1) << values[0].h << "   " << std::setprecision(4) << J(values[0].h, values[0].w, values[0].d) << "   " << values[0].w << "   " << values[0].d << "\n" << "\n";
    F_r3 = filter(F_noise, almin, r3);
    for (auto i : F_r3)
        std::cout << i << "\n";
    std::cout << "\n";


    values.clear();
    almin.clear();
    minD = 0;
    minW = 0;
    std::vector<double> F_r5;
    std::cout << " h    dist             alpha              w        d" << "\n";
    for (auto l = 0; l <= L; l++)
    {
        for (auto i = 0; i < N; i++)
        {
            al = alpha(r5);
            F_r5 = filter(F_noise, al, r5);
            W = w(F_r5);
            D = d(F_r5, F_noise);
            Dist = dist(W, D);

            if (Dist < minDist)
            {
                minDist = Dist;
                minW = W;
                minD = D;
                almin = al;
            }
        }
        H = (double)l / L;
        values.push_back({ H, minDist, almin, minW, minD });
        std::cout << std::fixed << std::setprecision(1) << H << "  " << std::setprecision(4) << minDist << "   [ ";
        for (auto i : almin)
            std::cout << i << " ";
        std::cout << "]  " << minW << "   " << minD << "\n";
        minDist = 100000;
    }

    std::sort(values.begin(), values.end(), [](auto x, auto y) {return x.dist < y.dist; });
    std::cout << "\n" << " h*     J        w        d" << "\n";
    std::cout << std::setprecision(1) << values[0].h << "   " << std::setprecision(4) << J(values[0].h, values[0].w, values[0].d) << "   " << values[0].w << "   " << values[0].d << "\n" << "\n";
    F_r5 = filter(F_noise, almin, r5);
    for (auto i : F_r5)
        std::cout << i << "\n";
    std::cout << "\n";
    
    return 0;
}

