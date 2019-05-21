#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

struct GlobalData
{
	double S;
	double K;
	double L;
	double alfa;      //Wsp�czynnik konwekcyjnej wymiany ciep�a
	double temp_ot;       //Temperatura otoczenia
	int nh;               //Liczba w�z��w
	double q;         //Pr�dko�� przekazywania ciep�a
	
	void ReadFile()
	{
		ifstream in("Data.txt");    //Wczytywanie z pliku
		in >> S >> K >> L >> alfa >> temp_ot >>nh >> q;
		in.close();
	}
};
GlobalData Data;

struct Element
{
	int ID[2];              //ID node'�w tworz�cych element
	double S;               //Pole przekroju
	double K;               //Wsp�czynnik przewodzenia ciep�a
	double L;               //D�ugo�� pr�ta
	double q;               //Pr�dko�� przekazywania ciep�a
	double HL[2][2];        //Macierz lokalna
	double PL;              //Wektor obci��e� lokalnych
};

struct Grid
{
	Element *Elements = new Element[Data.nh - 1];
};

void Elements(Grid grid)        //operacje na ka�dym z element�w
{
	for (int i = 0; i < Data.nh - 1; i++)
	{
		grid.Elements[i].ID[0] = i;
		grid.Elements[i].ID[1] = i + 1;
		grid.Elements[i].S = Data.S;
		grid.Elements[i].K = Data.K;
		grid.Elements[i].L = Data.L / (Data.nh - 1);
		grid.Elements[i].q = Data.q;

		for (int j = 0; j < 1; j++)
		{
			if (i < Data.nh - 2)
			{
				grid.Elements[i].HL[j][j] = (grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L;
				grid.Elements[i].HL[j][j + 1] = -(grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L;      //wzory w macierzy
				grid.Elements[i].HL[j + 1][j] = -(grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L;
				grid.Elements[i].HL[j + 1][j + 1] = (grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L;
			}
			else
			{
				grid.Elements[i].HL[j][j] = (grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L;
				grid.Elements[i].HL[j][j + 1] = -(grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L;
				grid.Elements[i].HL[j + 1][j] = -(grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L;
				grid.Elements[i].HL[j + 1][j + 1] = ((grid.Elements[i].S * grid.Elements[i].K) / grid.Elements[i].L) + (Data.alfa * grid.Elements[i].S);
			}
		}
	}
}

void Agregation(Grid grid)
{
	double **HG = new double *[Data.nh];
	for (int i = 0; i < Data.nh; i++)
	{
		HG[i] = new double[Data.nh];
		for (int j = 0; j < Data.nh; j++)
		{
			HG[i][j] = 0.0;
		}
	}

	double *PG = new double[Data.nh];
	for (int i = 0; i < Data.nh; i++)
	{
		PG[i] = 0.0;
	}

	double *TG = new double[Data.nh];
	for (int i = 0; i < Data.nh; i++)
	{
		TG[i] = 0.0;
	}

	for (int i = 0; i < Data.nh - 1; i++)
	{
		Element *temp = &grid.Elements[i];

		HG[temp->ID[0]][temp->ID[0]] = temp->HL[0][0];
		HG[temp->ID[0]][temp->ID[1]] = temp->HL[0][1];
		HG[temp->ID[1]][temp->ID[0]] = temp->HL[1][0];
		HG[temp->ID[1]][temp->ID[1]] = temp->HL[1][1];

		if (i>0)
		{
			HG[i][i] += HG[i][i];
		}
	}

	cout << "-- MACIERZ GLOBALNA HG --" << endl;;

	for (int i = 0; i < Data.nh; i++)
	{
		for (int j = 0; j < Data.nh; j++)
		{
			cout << HG[i][j] << "\t";
		}
		cout << " " << endl;
	}

	cout << endl;
	cout << "-- WEKTOR PG --" << endl;

	for (int i = 0; i <= Data.nh - 1; i++)
	{
		if (i == 0)
		{
			PG[i] = Data.q * Data.S;
		}
		else if (i == (Data.nh - 1))
		{
			PG[i] = -(Data.alfa * Data.temp_ot * Data.S);
		}
		else
		{
			PG[i] = 0;
		}
	}

	for (int i = 0; i < Data.nh; i++)
	{
		cout << PG[i];
		cout << endl;
	}

//Rozwi�zanie uk�adu r�wna�

	for (int k = 0; k<Data.nh - 1; k++)
	{
		for (int i = k + 1; i<Data.nh; i++)
			HG[i][k] /= HG[k][k];

		for (int i = k + 1; i<Data.nh; i++)
		for (int j = k + 1; j<Data.nh; j++)
			HG[i][j] -= HG[i][k] * HG[k][j];
	}

	double s;
	TG[0] = PG[0];

	for (int i = 1; i<Data.nh; i++)
	{
		s = 0;
		for (int j = 0; j<i; j++)
		s += HG[i][j] * TG[j];

		TG[i] = PG[i] - s;
	}

	TG[Data.nh - 1] /= HG[Data.nh - 1][Data.nh - 1];

	for (int i = Data.nh - 2; i >= 0; i--)
	{
		s = 0;
		for (int j = i + 1; j<Data.nh; j++) s += HG[i][j] * TG[j];
		TG[i] = (TG[i] - s) / HG[i][i];
	}

	cout << endl;
	cout << "** WYNIKI **" << endl;
	for (int i = 0; i<Data.nh; i++)
	{
		cout << "T" << i+1 << " = " << TG[i] * -1 << endl;
	}
}

int main()
{
	Data.ReadFile();
	Grid grid;
	Elements(grid);

	Agregation(grid);
	
	cout << endl;
	system("pause");
}