//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Unit1.h"
#include <cmath>
#include <vector>
//#include "shared.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "VCLTee.TeeSurfa"
#pragma link "VCLTee.TeeSurfa"
#pragma resource "*.dfm"
TForm1 *Form1;
bool stop=true;
const double _dx = .05, _dy = .05, _dt = _dx / 10.; //, _dt;

	size_t N_x = 300;
	size_t N_y = 60;
	std::vector<std::vector<double> > psi(N_x, std::vector<double>(N_y));
	std::vector<std::vector<double> > tet(N_x, std::vector<double>(N_y));
	std::vector<std::vector<double> > e(N_x, std::vector<double>(N_y));
	std::vector<std::vector<double> > u(N_x, std::vector<double>(N_y));
	std::vector<std::vector<double> > v(N_x, std::vector<double>(N_y));



	double Re = 100.;
	double Ri = 10.;
	const double Pe = Re;

	double w = 2. / (1 + sin(M_PI * (_dx + _dy) / 2.));

    const size_t iter_MAX = 5000;
    const double eps = 1e-2;

    size_t iter_count = 0;
    size_t iter_best = iter_MAX;
    double E = 0;
    double q = pow(_dx / _dy, 2);
	double w_best = 0;
	double w_step = 15e-3;
    size_t new_tet = 0;
	 // начальные
        const double _p0 = 10., _x0 = _dx * N_x / 2., _r0 = .06;
	   //	shared::randomize();


	size_t ct = 1;
	void start()
	{
	   for (size_t i = 0; i < N_x; i++)
	   {
		for (size_t j = 0; j < N_y; j++)
		{
			psi[i][j] = 0;
			tet[i][j] = 0;
			e[i][j] = 0;
			u[i][j] = 0;
			v[i][j] = 0;
		}
	   }
	   for (size_t i = 0; i < N_x; i++)
	   {
			tet[i][0] = (double)rand() / RAND_MAX;
		}
    }
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void __fastcall TForm1::StartClick(TObject *Sender)
{
if (stop)
{
   stop = false;
   Start->Caption = "Stop";
}
else
{
   stop = true;
   Start->Caption = "Continue";
}
	while (!stop) {
    ct++;
        { // адвекция диффузия схема против потока  e
            for (size_t i = 1; i < N_x - 1; i++) {
                for (size_t j = 1; j < N_y - 1; j++) {
                    unsigned a = 1, b = 1;
                    if (u[i][j] > 0) {
						a = 0;
                    }
                    if (v[i][j] > 0) {
                        b = 0;
                    }
                    e[i][j] =
						e[i][j] -
                        _dt / _dx * u[i][j] * (1 - a) *
                            (e[i][j] - e[i - 1][j]) -
                        _dt / _dx * a * u[i][j] * (e[i + 1][j] - e[i][j]) -
                        _dt / _dy * v[i][j] * (1 - b) *
                            (e[i][j] - e[i][j - 1]) -
						_dt / _dy * b * v[i][j] * (e[i][j + 1] - e[i][j]) +
                        _dt / Re *
                            ((e[i + 1][j] - 2. * e[i][j] + e[i - 1][j]) / _dx /
                                    _dx +
                                (e[i][j + 1] - 2. * e[i][j] + e[i][j - 1]) /
                                    _dy / _dy) -
						Ri * _dt / 2. / _dx * (tet[i + 1][j] - tet[i - 1][j]);
                }
            }
            // края
            for (size_t j = 1; j < N_y - 1; j++) {
                unsigned a = 1, b = 1;
				if (u[0][j] > 0) {
                    a = 0;
                }
                if (v[0][j] > 0) {
                    b = 0;
                }
				e[0][j] =
                    e[0][j] -
                    _dt / _dx * u[0][j] * (1 - a) * (e[0][j] - e[N_x - 2][j]) -
                    _dt / _dx * a * u[0][j] * (e[0 + 1][j] - e[0][j]) -
                    _dt / _dy * v[0][j] * (1 - b) * (e[0][j] - e[0][j - 1]) -
                    _dt / _dy * b * v[0][j] * (e[0][j + 1] - e[0][j]) +
					_dt / Re *
                        ((e[0 + 1][j] - 2. * e[0][j] + e[N_x - 2][j]) / _dx /
                                _dx +
                            (e[0][j + 1] - 2. * e[0][j] + e[0][j - 1]) / _dy /
                                _dy) -
                    Ri * _dt / 2. / _dx * (tet[0 + 1][j] - tet[N_x - 2][j]);
			}

            for (size_t j = 1; j < N_y - 1; j++) {
                e[N_x - 1][j] = e[0][j];
            }
		}

        { // адвекция диффузия схема против потока tet
            for (size_t i = 1; i < N_x - 1; i++) {
                for (size_t j = 1; j < N_y - 1; j++) {
                    unsigned a = 1, b = 1;
					if (u[i][j] > 0) {
                        a = 0;
                    }
                    if (v[i][j] > 0) {
                        b = 0;
                    }
					tet[i][j] =
                        tet[i][j] -
                        _dt / _dx * u[i][j] * (1 - a) *
                            (tet[i][j] - tet[i - 1][j]) -
                        _dt / _dx * a * u[i][j] * (tet[i + 1][j] - tet[i][j]) -
                        _dt / _dy * v[i][j] * (1 - b) *
							(tet[i][j] - tet[i][j - 1]) -
                        _dt / _dy * b * v[i][j] * (tet[i][j + 1] - tet[i][j]) +
                        _dt / Pe *
                            ((tet[i + 1][j] - 2. * tet[i][j] + tet[i - 1][j]) /
                                    _dx / _dx +
                                (tet[i][j + 1] - 2. * tet[i][j] +
									tet[i][j - 1]) /
                                    _dy / _dy);
                }
            }
            // края
            for (size_t j = 1; j < N_y - 1; j++) {
				unsigned a = 1, b = 1;
                if (u[0][j] > 0) {
                    a = 0;
                }
                if (v[0][j] > 0) {
                    b = 0;
				}
                tet[0][j] =
                    tet[0][j] -
                    _dt / _dx * u[0][j] * (1 - a) *
                        (tet[0][j] - tet[N_x - 2][j]) -
                    _dt / _dx * a * u[0][j] * (tet[0 + 1][j] - tet[0][j]) -
					_dt / _dy * v[0][j] * (1 - b) *
                        (tet[0][j] - tet[0][j - 1]) -
                    _dt / _dy * b * v[0][j] * (tet[0][j + 1] - tet[0][j]) +
                    _dt / Pe *
                        ((tet[0 + 1][j] - 2. * tet[0][j] + tet[N_x - 2][j]) /
                                _dx / _dx +
							(tet[0][j + 1] - 2. * tet[0][j] + tet[0][j - 1]) /
                                _dy / _dy);
            }

            for (size_t j = 1; j < N_y - 1; j++) {
                tet[N_x - 1][j] = tet[0][j];
			}
        }

        { //элиптическое
            iter_count = 0;
            do {
				E = 0;
                for (size_t i = 1; i < N_x - 1; i++) {
                    for (size_t j = 1; j < N_y - 1; j++) {
                        double prev = psi[i][j];
                        psi[i][j] = (psi[i + 1][j] + psi[i - 1][j] +
                                        q * (psi[i][j + 1] + psi[i][j - 1])) /
										((double)(2. + 2. * q)) -
                                    _dx * _dx / (2. + 2. * q) * e[i][j];
                        psi[i][j] = prev + w * (psi[i][j] - prev);
                        double er = abs(prev - psi[i][j]);
                        if (E < er) {
                            E = er;
						}
                    }
                }
                iter_count++;
            } while (eps < E && iter_count < iter_best);
            if (iter_count < iter_best) {
				iter_best = iter_count;
            } else {
                iter_best = iter_MAX;
                w += w_step;
                if (w > 1.9) {
                    w = 1;
				}
            }
            // для краев
            size_t iter_count2 = 0;
            do {
                E = 0;
				for (size_t j = 1; j < N_y - 1; j++) {
                    double prev = psi[0][j];
                    psi[0][j] = (psi[1][j] + psi[N_y - 2][j] +
                                    q * (psi[0][j + 1] + psi[0][j - 1])) /
                                    ((double)(2. + 2. * q)) -
                                _dx * _dx / (2. + 2. * q) * e[0][j];
					psi[0][j] = prev + w * (psi[0][j] - prev);
                    double er = abs(prev - psi[0][j]);
                    if (E < er) {
                        E = er;
                    }
                }
				iter_count2++;
            } while (eps < E && iter_count2 < iter_MAX);
            for (size_t j = 1; j < N_y - 1; j++) {
                psi[N_x - 1][j] = psi[0][j];
            }
        }

        { // пересчет поля скоростей
            for (size_t i = 0; i < N_x; i++) {
                for (size_t j = 1; j < N_y - 1; j++) {
                    u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / 2. / _dy;
                }
			}
            for (size_t i = 1; i < N_x - 1; i++) {
                for (size_t j = 0; j < N_y; j++) {
                    v[i][j] = -1 * (psi[i + 1][j] - psi[i - 1][j]) / 2. / _dx;
                }
            }

            for (size_t j = 0; j < N_y; j++) {
                v[0][j] = -1 * (psi[1][j] - psi[N_x - 2][j]) / 2. / _dx;
                v[N_x - 1][j] = v[0][j];
            }
        }

        { //                    пересчет
            { // стена, твердая граница
                for (size_t i = 0; i < N_x; i++) {
                    tet[i][N_y - 1] = tet[i][N_y - 2] * 0.7;
                    u[i][N_y - 1] = 0;
					u[i][0] = 0;
                    v[i][N_y - 1] = 0;
                    v[i][0] = 0;
                    psi[i][N_y - 1] = 0;
                    psi[i][0] = 0;
                    e[i][N_y - 1] = 2. * (psi[i][N_y - 2]) / _dy / _dy;
					e[i][0] = 2. * (psi[i][1]) / _dy / _dy;
                }
            };
        }
		{ // вывод

				Form1->Series1->Clear();
                for (size_t i = 0; i < N_x; i++) {
                    for (size_t j = 0; j < N_y; j++) {
                        Form1->Series1->AddXYZ(i, tet[i][j], j);
                    }
                }
				Form1->Caption =
                    "Re" + FloatToStr(Re) + " | Ri" + FloatToStr(Ri) + " | Pe" +
                    FloatToStr(Pe) + " | t:" + FloatToStr(ct);
				Form1->Chart1->Repaint();
            Application->ProcessMessages();
        }
	}
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button2Click(TObject *Sender)
{
   ct=1;
   Re = StrToFloat(LabeledEdit1->Text);
   Ri = StrToFloat(LabeledEdit2->Text);
   start();
   stop = true;
   Start->Caption = "Start";
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

