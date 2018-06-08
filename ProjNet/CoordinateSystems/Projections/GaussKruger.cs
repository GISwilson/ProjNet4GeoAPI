using System;
using System.Collections.Generic;
using GeoAPI.CoordinateSystems;
using GeoAPI.CoordinateSystems.Transformations;

namespace ProjNet.CoordinateSystems.Projections
{
    internal class GaussKruger:MapProjection
    {
        //private bool is3dZone;// is 3 degree Gauss Kruger
        //private bool containZone;//is contain zone num in x coordinate
        //private int zoneNum;

        //private double f;//曲率
        private double e1;//第一离心率
        private double e2;//第二离心率
        private double n;//卯酉曲率半径
        //private const double PiPerDegree = Math.PI / 180;

        public GaussKruger(IEnumerable<ProjectionParameter> parameters, MapProjection inverse) : base(parameters, inverse)
        {
            var a = _semiMajor;
            var b = _semiMinor;
            e1 = (a * a - b * b) / (a * a);
            e2 = (a * a - b * b) / (b * b);
            n = (1.0 - Math.Sqrt(1.0 - this.e1)) / (1.0 + Math.Sqrt(1.0 - this.e1));
            Authority = "EPSG";
            //is3dZone = ((int) central_meridian) % 3 != 0;
            //containZone = false_easting - 500000 > 1e-5;
            //if (containZone)
            //{
            //    zoneNum = Convert.ToInt16(Math.Round(false_easting - 500000) / 1000000.0);
            //}
        }

        public GaussKruger(IEnumerable<ProjectionParameter> parameters) : base(parameters)
        {
            var a = _semiMajor;
            var b = _semiMinor;
            e1 = (a * a - b * b) / (a * a);
            e2 = (a * a - b * b) / (b * b);
            n = (1.0 - Math.Sqrt(1.0 - this.e1)) / (1.0 + Math.Sqrt(1.0 - this.e1));
            Authority = "EPSG";
            //is3dZone = ((int) central_meridian) % 3 != 0;
            //containZone = false_easting - 500000 > 1e-5;
            //if (containZone)
            //{
            //    zoneNum = Convert.ToInt16(Math.Round(false_easting - 500000) / 1000000.0);
            //}
        }

        public override IMathTransform Inverse()
        {
            if (_inverse == null)
                _inverse = new GaussKruger(_Parameters.ToProjectionParameter(), this);
            return _inverse;
        }

        //正算
        protected override double[] RadiansToMeters(double[] lonlat)
        {
            if (double.IsNaN(lonlat[0]) || double.IsNaN(lonlat[1]))
                return new[] { Double.NaN, Double.NaN, };

            var longitude0 = central_meridian;//换算为弧度的中央经线
            var longitude = lonlat[0];
            var latitude = lonlat[1];
            var longitude1 = longitude - longitude0; //经度转换为弧度
            var latitude1 = latitude;//不必换算了latitude * PiPerDegree; //纬度转换为弧度
            var n = _semiMajor / Math.Sqrt(1.0 - e1 * Math.Sin(latitude1) * Math.Sin(latitude1));
            var tan = Math.Tan(latitude1);
            var tao = Math.Sqrt(e2) * Math.Cos(latitude1);
            var cos = Math.Cos(latitude1);
            var sin = Math.Sin(latitude1);
            var cos2 = cos * cos;

            var a0 = _semiMajor * ((1 - e1 / 4 - 3 * e1 * e1 / 64 - 5 * e1 * e1 * e1 / 256) * latitude1 - (3 * e1 / 8 + 3 * e1 * e1 / 32 + 45 * e1 * e1
                                                                                                          * e1 / 1024) * Math.Sin(2 * latitude1)
                                  + (15 * e1 * e1 / 256 + 45 * e1 * e1 * e1 / 1024) * Math.Sin(4 * latitude1) - (35 * e1 * e1 * e1 / 3072) * Math.Sin(6 * latitude1));
            var a1 = n * cos;
            var a2 = n * tan * cos2 / 2;
            var a3 = n * (1 - tan * tan + tao * tao) * cos2 * cos / 6;
            var a4 = n * tan * (5 - tan * tan + 9 * tao * tao + 4 * tao * tao * tao * tao) * cos2 * cos2 / 24;
            //var a5 = NN * (5 - 18 * T * T + T * T * T * T + 14 * tao * tao - 58 * T * T * tao * tao) * cos2 * cos2 *
            //         cos / 120;//simple method
            var a5 = n * cos / 120 * (1 - 20 * cos2 + 24 * cos2 * cos2 + (-58 * cos2 + 72 * cos2 * cos2) * tao * tao +
                                   (-64 * cos2 + 77 * cos2 * cos2) * tao * tao * tao * tao +
                                   (-24 * cos2 + 28 * cos2 * cos2) * tao * tao * tao * tao * tao * tao);
            //var a6 = NN * T * (61 - 58 * T * T + T * T * T * T + 270 * tao * tao - 330 * T * T * tao * tao) * cos2 *
            //         cos2 * cos2 / 720;//simple method
            var a6 = n * sin * cos / 720 * (1 - 60 * cos2 + 120 * cos2 * cos2 +
                                         (-330 * cos2 + 600 * cos2 * cos2) * tao * tao +
                                         (-680 * cos2 + 1125 * cos2 * cos2) * tao * tao * tao * tao +
                                         (-600 * cos2 + 924 * cos2 * cos2) * tao * tao * tao * tao * tao * tao +
                                         (-192 * cos2 + 280 * cos2 * cos2) * tao * tao * tao * tao * tao * tao * tao *
                                         tao);
            var a7 = n * cos / 5040 * (-1 + 182 * cos2 - 840 * cos2 * cos2 + 720 * cos2 * cos2 * cos2 +
                                        (1771 * cos2 * cos2 - 6840 * cos2 * cos2 * cos2 +
                                         5400 * cos2 * cos2 * cos2 * cos2) * tao * tao);
            

            var yval = a0 + a2 * longitude1 * longitude1 + a4 * longitude1 * longitude1 * longitude1 * longitude1 +
                       a6 * longitude1 * longitude1 * longitude1 * longitude1 * longitude1 * longitude1;
            var xval = a1 * longitude1 + a3 * longitude1 * longitude1 * longitude1 +
                       a5 * longitude1 * longitude1 * longitude1 * longitude1 * longitude1 +
                       a7 * longitude1 * longitude1 * longitude1 * longitude1 * longitude1 * longitude1 * longitude1;

            //var x0 = 0d;
            //if (containZone)
            //{
            //    x0 = 1000000L * zoneNum + 500000L;
            //}
            //xval = xval + x0;
            return new double[] { xval, yval };
        }

        //反算
        protected override double[] MetersToRadians(double[] p)
        {
            /* Inverse equations
			  -----------------*/
            double x = p[0]; //* _metersPerUnit - this._falseEasting;
            double y = p[1];// * _metersPerUnit - this._falseNorthing;

            
            
            double num12 = y / (_semiMajor * (((1.0 - (e1 / 4.0)) - (((3.0 * e1) * e1) / 64.0)) - ((((5.0 * e1) * e1) * e1) / 256.0)));
            double d = (((num12 + ((((3.0 * n) / 2.0) - ((((27.0 * n) * n) * n) / 32.0)) * Math.Sin(2.0 * num12))) + (((((21.0 * n) * n) / 16.0) - (((((55.0 * n) * n) * n) * n) / 32.0)) * Math.Sin(4.0 * num12))) + (((((151.0 * n) * n) * n) / 96.0) * Math.Sin(6.0 * num12))) + ((((((1097.0 * n) * n) * n) * n) / 512.0) * Math.Sin(8.0 * num12));
            double num14 = (e2 * Math.Cos(d)) * Math.Cos(d);
            double num15 = Math.Tan(d) * Math.Tan(d);
            double num16 = _semiMajor / Math.Sqrt(1.0 - ((e1 * Math.Sin(d)) * Math.Sin(d)));
            double num17 = (_semiMajor * (1.0 - e1)) / Math.Sqrt(((1.0 - ((e1 * Math.Sin(d)) * Math.Sin(d))) * (1.0 - ((e1 * Math.Sin(d)) * Math.Sin(d)))) * (1.0 - ((e1 * Math.Sin(d)) * Math.Sin(d))));
            double num18 = x / num16;
            double num19 = central_meridian + (((num18 - ((((((1.0 + (2.0 * num15)) + num14) * num18) * num18) * num18) / 6.0)) + (((((((((((5.0 - (2.0 * num14)) + (28.0 * num15)) - ((3.0 * num14) * num14)) + (8.0 * e2)) + ((24.0 * num15) * num15)) * num18) * num18) * num18) * num18) * num18) / 120.0)) / Math.Cos(d));
            double num20 = d - (((num16 * Math.Tan(d)) / num17) * ((((num18 * num18) / 2.0) - (((((((((5.0 + (3.0 * num15)) + (10.0 * num14)) - ((4.0 * num14) * num14)) - (9.0 * e2)) * num18) * num18) * num18) * num18) / 24.0)) + ((((((((((((61.0 + (90.0 * num15)) + (298.0 * num14)) + ((45.0 * num15) * num15)) - (256.0 * e2)) - ((3.0 * num14) * num14)) * num18) * num18) * num18) * num18) * num18) * num18) / 720.0)));
            //var longitude = num19 / PiPerDegree;//不必换算
            //var latitude = num20 / PiPerDegree;
            return new double[] { num19, num20 };
        }
    }
}
