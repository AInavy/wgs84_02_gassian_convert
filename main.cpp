#include <iostream>
#include <cmath>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>



using namespace std;
#define PAI 3.1415926


namespace strtool
{
    string trim(const string& str)
    {
        string::size_type pos = str.find_first_not_of(' ');
        if (pos == string::npos)
        {
            return str;
        }
        string::size_type pos2 = str.find_last_not_of(' ');
        if (pos2 != string::npos)
        {
            return str.substr(pos, pos2 - pos + 1);
        }
        return str.substr(pos);
    }
    int split(const string& str, vector<string>& ret_, string sep = ",")
    {
        if (str.empty())
        {
            return 0;
        }
        string tmp;
        string::size_type pos_begin = str.find_first_not_of(sep);
        string::size_type comma_pos = 0;
        while (pos_begin != string::npos)
        {
            comma_pos = str.find(sep, pos_begin);
            if (comma_pos != string::npos)
            {
                tmp = str.substr(pos_begin, comma_pos - pos_begin);
                pos_begin = comma_pos + sep.length();
            }
            else
            {
                tmp = str.substr(pos_begin);
                pos_begin = comma_pos;
            }
            if (!tmp.empty())
            {
                ret_.push_back(tmp);
                tmp.clear();
            }
        }
        return 0;
    }
};

//#define element_of(o) (sizeof(o) / sizeof(o[0]))
//#define UNUSED(o) (o = o)

typedef struct iPoint
{
    int x,y;

    iPoint()
    {
        x = y = 0;
    }
} PointI;

typedef struct dPoint
{
    double x,y;
} PointD;

typedef struct dPoint3d
{
    double x,y,z;
} Point3D;



double a = 6378245.0;
double ee = 0.00669342162296594323;
double transformLat(double x, double y)
{
    double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
    ret += (20.0 * sin(6.0 * x * PAI) + 20.0 * sin(2.0 * x * PAI)) * 2.0 / 3.0;
    ret += (20.0 * sin(y * PAI) + 40.0 * sin(y / 3.0 * PAI)) * 2.0 / 3.0;
    ret += (160.0 * sin(y / 12.0 * PAI) + 320 * sin(y * PAI / 30.0)) * 2.0 / 3.0;
    return ret;
}

double transformLon(double x, double y)
{
    double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
    ret += (20.0 * sin(6.0 * x * PAI) + 20.0 * sin(2.0 * x * PAI)) * 2.0 / 3.0;
    ret += (20.0 * sin(x * PAI) + 40.0 * sin(x / 3.0 * PAI)) * 2.0 / 3.0;
    ret += (150.0 * sin(x / 12.0 * PAI) + 300.0 * sin(x / 30.0 * PAI)) * 2.0 / 3.0;
    return ret;
}

void transform(double lat, double lon, double& _relat, double& _relon)
{
    double dLat = transformLat(lon - 105.0, lat - 35.0);
    double dLon = transformLon(lon - 105.0, lat - 35.0);
    double radLat = lat / 180.0 * PAI;
    double magic = sin(radLat);
    magic = 1 - ee * magic * magic;
    double sqrtMagic = sqrt(magic);
    dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * PAI);
    dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * PAI);
    _relat = lat + dLat;
    _relon = lon + dLon;
}

void gps84_To_Gcj02(dPoint& wgs84, dPoint& gcj02)
{
    double lat = wgs84.y;
    double lon = wgs84.x;
    double latout, lonout;
    double dLat = transformLat(lon - 105.0, lat - 35.0);
    double dLon = transformLon(lon - 105.0, lat - 35.0);
    double radLat = lat / 180.0 * PAI;
    double magic = sin(radLat);
    magic = 1 - ee * magic * magic;
    double sqrtMagic = sqrt(magic);
    dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * PAI);
    dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * PAI);
    latout = lat + dLat;
    lonout = lon + dLon;
    gcj02.x = lonout;
    gcj02.y = latout;
}

void gcj02_To_Gps84(dPoint& _va02, dPoint& _va84)
{
    double lat = _va02.y;
    double lon = _va02.x;
    double retlat, retlon;
    transform(lat, lon, retlat, retlon);
    double lontitude = lon * 2 - retlon;
    double latitude = lat * 2 - retlat;
    _va84.x = lontitude;
    _va84.y = latitude;
}

dPoint GaussProjInvCal(dPoint xoy)
{

    double X = xoy.x;
    double Y = xoy.y;

    int ProjNo; int ZoneWide; ////?????
    double longitude1, latitude1, longitude0, X0, Y0, xval, yval;
    double e1, e2, f, a, ee, NN, T, C, M, D, R, u, fai, d_PAI;
    d_PAI = 0.0174532925199433; ////3.1415926535898/180.0;
    a = 6378245.0; f = 1.0 / 298.3; //54??????????????????????
    ////a=6378140.0; f=1/298.257; //80?????????????????????
    ZoneWide = 6; ////6???????
    ProjNo = (int)(X / 1000000L); //?????????
    longitude0 = (ProjNo - 1) * ZoneWide + ZoneWide / 2;
    longitude0 = longitude0 * d_PAI; //????????
    X0 = ProjNo * 1000000L + 500000L;
    Y0 = 0;
    xval = X - X0; yval = Y - Y0; //????????????????
    e2 = 2 * f - f*f;
    e1 = (1.0 - sqrt(1 - e2)) / (1.0 + sqrt(1 - e2));
    ee = e2 / (1 - e2);
    M = yval;
    u = M / (a*(1 - e2 / 4 - 3 * e2*e2 / 64 - 5 * e2*e2*e2 / 256));
    fai = u + (3 * e1 / 2 - 27 * e1*e1*e1 / 32)*sin(2 * u) + (21 * e1*e1 / 16 - 55 * e1*e1*e1*e1 / 32)*sin(
            4 * u)
          + (151 * e1*e1*e1 / 96)*sin(6 * u) + (1097 * e1*e1*e1*e1 / 512)*sin(8 * u);
    C = ee*cos(fai)*cos(fai);
    T = tan(fai)*tan(fai);
    NN = a / sqrt(1.0 - e2*sin(fai)*sin(fai));
    R = a*(1 - e2) / sqrt((1 - e2*sin(fai)*sin(fai))*(1 - e2*sin(fai)*sin(fai))*(1 - e2*sin(fai)*sin(fai)));
    D = xval / NN;
    //????????(Longitude) ????(Latitude)
    longitude1 = longitude0 + (D - (1 + 2 * T + C)*D*D*D / 6 + (5 - 2 * C + 28 * T - 3 * C*C + 8 * ee + 24 * T*T)*D
                                                               *D*D*D*D / 120) / cos(fai);
    latitude1 = fai - (NN*tan(fai) / R)*(D*D / 2 - (5 + 3 * T + 10 * C - 4 * C*C - 9 * ee)*D*D*D*D / 24
                                         + (61 + 90 * T + 298 * C + 45 * T*T - 256 * ee - 3 * C*C)*D*D*D*D*D*D / 720);

    // 	//????????? DD
    // 	*longitude = longitude1 / iPAI;
    // 	*latitude = latitude1 / iPAI;

    dPoint bol;
    bol.x = longitude1 / d_PAI;
    bol.y = latitude1 / d_PAI;

    return bol;
}


dPoint GaussProjCal(dPoint bol)
{
    double longitude = bol.x;//????
    double latitude = bol.y; //????

    int ProjNo = 0; int ZoneWide; ////?????
    double longitude1, latitude1, longitude0, latitude0, X0, Y0, xval, yval;
    double a, f, e2, ee, NN, T, C, A, M, d_PAI;
    d_PAI = 0.0174532925199433; ////3.1415926535898/180.0;
    ZoneWide = 6; ////6???????
    a = 6378245.0; f = 1.0 / 298.3; //54??????????????????????
    ////a=6378140.0; f=1/298.257; //80?????????????????????
    ProjNo = (int)(longitude / ZoneWide);
    longitude0 = ProjNo * ZoneWide + ZoneWide / 2;
    longitude0 = longitude0 * d_PAI;
    latitude0 = 0;
    longitude1 = longitude * d_PAI; //???????????????
    latitude1 = latitude * d_PAI; //???????????????
    e2 = 2 * f - f*f;
    ee = e2*(1.0 - e2);
    NN = a / sqrt(1.0 - e2*sin(latitude1)*sin(latitude1));
    T = tan(latitude1)*tan(latitude1);
    C = ee*cos(latitude1)*cos(latitude1);
    A = (longitude1 - longitude0)*cos(latitude1);

    M = a*((1 - e2 / 4 - 3 * e2*e2 / 64 - 5 * e2*e2*e2 / 256)*latitude1 - (3 * e2 / 8 + 3 * e2*e2 / 32 + 45 * e2*e2*e2 / 1024)*sin(2 * latitude1)
           + (15 * e2*e2 / 256 + 45 * e2*e2*e2 / 1024)*sin(4 * latitude1) - (35 * e2*e2*e2 / 3072)*sin(6 * latitude1));
    xval = NN*(A + (1 - T + C)*A*A*A / 6 + (5 - 18 * T + T*T + 72 * C - 58 * ee)*A*A*A*A*A / 120);
    yval = M + NN*tan(latitude1)*(A*A / 2 + (5 - T + 9 * C + 4 * C*C)*A*A*A*A / 24
                                  + (61 - 58 * T + T*T + 600 * C - 330 * ee)*A*A*A*A*A*A / 720);
    X0 = 1000000L * (ProjNo + 1) + 500000L;
    Y0 = 0;
    xval = xval + X0; yval = yval + Y0;

    // 	*X = xval;
    // 	*Y = yval;

    dPoint xoy;
    xoy.x = xval;
    xoy.y = yval;

    return xoy;

}

double getDistance(const dPoint& p1, const dPoint& p2)
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;

    return sqrt(dx*dx + dy* dy);
}



int main()
{
    ifstream fin("/home/mahj/PycharmProjects/json_parse/duantoulu_obj_filter.txt");
    ofstream fout("/home/mahj/PycharmProjects/json_parse/duantoulu_obj_filter_WGS84.txt");
    string s;
    while(getline(fin,s))
    {
        vector<string> s_vec;
        strtool::split(s,s_vec," ");
        double G = atof(s_vec[1].c_str());
        double I = atof(s_vec[2].c_str());

        dPoint had02;
        had02 = {G,I};

        dPoint had84;
        gcj02_To_Gps84(had02,had84);



        //cout<<s_vec[0]<<" ";
        cout<<setiosflags(ios::fixed)<<setprecision(8)<<had84.x<<" "<<had84.y<<" "<<s_vec[3]<<endl;

        //fout<<s_vec[0]<<" ";
        fout<<setiosflags(ios::fixed)<<setprecision(8)<<had84.x<<" "<<had84.y<<" "<<s_vec[3]<<endl;

        s_vec.clear();
    }



    return 0;
}
