//////////////////////////////////////////////////////////////////////////
///////  系统信息的数据结构    /////////////////////////////
//////////////////////////////////////////////////////////////////////////

#ifndef  _GENERALINFO_H
#define _GENERALINFO_H
 struct SystemInfo    ///system information
{
	int SystemTotalBusNum;    ///Total bus number
	int GeneralLineNum;       ///normal line sequence number
	int LandBranchNum;        //land branch amount
	int TransformerNum;       //transformer amount
	double BasePower;         ///base power
	double BaseVoltage;      ///base voltage
	int SwingBusNo;          ///balance bus number
	int MultiLineNum;        ///multiline number
};
struct LineInfo   ///Line information
{
     int BusINo;     ///i bus number
	 int BusJNo;     ///j bus number
	 double R;       ///line resistance
	 double X;       //line impedance
	 double Yk;      //line ground impedance
};
struct TransformerInfo     ///transformer line information
{
     int BusINo;    ///i bus number
	 int BusJNo;    ////none-unit side number
     double R;      //transformer resistance
	 double X;      //transformer impedance
     double Ratio;  ///transformer ratio
};
struct BusInfo
{
      int BusNo;     ///bus number
	  int BusType;   ///bus type
	  double PG;      //generator active power output
	  double QG;      //generator reactive power output
	  double PD;      //active power demand
	  double QD;      //reactive power demand
      double Voltage_e;  ///real part of bus voltage
	  double Voltage_f;  ///imaginary part of bus voltage
	  double VolMod;     //voltage mode
};
#endif
