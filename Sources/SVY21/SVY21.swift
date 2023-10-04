// The Swift Programming Language
// https://docs.swift.org/swift-book

import Foundation

final class SVY21 {
    // WGS84 Datum
    static private let a:Double = 6378137
    static private let f:Double = 1.0 / 298.257223563
    static private let PI:Double = atan(1.0)*4.0
    
    // SVY21 Projection
    // Fundamental point: Base 7 at Pierce Resevoir.
    // Latitude: 1 22 02.9154 N, longitude: 103 49 31.9752 E (of Greenwich).
    
    // Known Issue: Setting (oLat, oLon) to the exact coordinates specified above
    // results in computation being slightly off. The values below give the most
    // accurate represenation of test data.
    static private let oLat:Double = 1.366666;     // origin's lat in degrees
    static private let oLon:Double = 103.833333;   // origin's lon in degrees
    static private let oN:Double = 38744.572;      // false Northing
    static private let oE:Double = 28001.642;      // false Easting
    static private let k:Double = 1.0;                // scale factor
    
    static private let b:Double = a * (1.0 - f);
    static private let e2:Double = (2.0 * f) - (f * f);
    static private let e4:Double = e2 * e2;
    static private let e6:Double = e4 * e2;
    static private let A0:Double = 1.0 - (e2 / 4.0) - (3.0 * e4 / 64.0) - (5.0 * e6 / 256.0);
    static private let A2:Double = (3.0 / 8.0) * (e2 + (e4 / 4.0) + (15.0 * e6 / 128.0));
    static private let A4:Double = (15.0 / 256.0) * (e4 + (3.0 * e6 / 4.0));
    static private let A6:Double = 35.0 * e6 / 3072.0;
    
    private static func calcM(_ lat: Double) -> Double {
        
        let latR:Double = lat * PI / 180.0;
        return a * ((A0 * latR) - (A2 * sin(2 * latR)) + (A4 * sin(4 * latR)) - (A6 * sin(6 * latR)));
    }
    
    private static func calcRho(_ sin2Lat: Double) -> Double {
        
        let num:Double = a * (1 - e2);
        let denom:Double = pow(1 - e2 * sin2Lat, 3.0 / 2.0);
        return num / denom;
    }
    
    private static func calcV(_ sin2Lat: Double) -> Double {
        
        let poly:Double = 1.0 - e2 * sin2Lat;
        return a / sqrt(poly);
    }
    
    //    LatLon to SVY21
    static func latLonToSVY21(_ lat: Double, _ lon: Double) -> (northing: Double, easting:Double) {
        //Returns a pair (N, E) representing Northings and Eastings in SVY21.
        
        let latR:Double = lat * PI / 180.0;
        let sinLat:Double = sin(latR);
        let sin2Lat:Double = sinLat * sinLat;
        let cosLat:Double = cos(latR);
        let cos2Lat:Double = cosLat * cosLat;
        let cos3Lat:Double = cos2Lat * cosLat;
        let cos4Lat:Double = cos3Lat * cosLat;
        let cos5Lat:Double = cos4Lat * cosLat;
        let cos6Lat:Double = cos5Lat * cosLat;
        let cos7Lat:Double = cos6Lat * cosLat;
        
        let rho:Double = calcRho(sin2Lat);
        let v:Double = calcV(sin2Lat);
        let psi:Double = v / rho;
        let t:Double = tan(latR);
        let w:Double = (lon - oLon) * PI / 180.0;
        
        let M:Double = calcM(lat);
        let Mo:Double = calcM(oLat);
        
        let w2:Double = w * w;
        let w4:Double = w2 * w2;
        let w6:Double = w4 * w2;
        let w8:Double = w6 * w2;
        
        let psi2:Double = psi * psi;
        let psi3:Double = psi2 * psi;
        let psi4:Double = psi3 * psi;
        
        let t2:Double = t * t;
        let t4:Double = t2 * t2;
        let t6:Double = t4 * t2;
        
        //    Compute Northing
        let nTerm1:Double = w2 / 2.0 * v * sinLat * cosLat;
        let nTerm2:Double = w4 / 24.0 * v * sinLat * cos3Lat * (4.0 * psi2 + psi - t2);
        let nTerm3:Double = w6 / 720.0 * v * sinLat * cos5Lat * ((8.0 * psi4) * (11.0 - 24.0 * t2) - (28.0 * psi3) * (1.0 - 6.0 * t2) + psi2 * (1.0 - 32.0 * t2) - psi * 2.0 * t2 + t4);
        let nTerm4:Double = w8 / 40320.0 * v * sinLat * cos7Lat * (1385.0 - 3111.0 * t2 + 543.0 * t4 - t6);
        
        let northing:Double = oN + k * (M - Mo + nTerm1 + nTerm2 + nTerm3 + nTerm4);
        
        //    Compute Easting
        let eTerm1:Double = w2 / 6.0 * cos2Lat * (psi - t2);
        let eTerm2:Double = w4 / 120.0 * cos4Lat * ((4.0 * psi3) * (1.0 - 6.0 * t2) + psi2 * (1.0 + 8.0 * t2) - psi * 2.0 * t2 + t4);
        let eTerm3:Double = w6 / 5040.0 * cos6Lat * (61.0 - 479.0 * t2 + 179.0 * t4 - t6);
        
        let easting:Double = oE + k * v * w * cosLat * (1.0 + eTerm1 + eTerm2 + eTerm3);
        
        return (northing, easting)
    }
    
    
    //    SVY21 to LatLon
    static func SVY21ToLatLon(_ northing: Double, _ easting: Double) -> (lat:Double, lon:Double) {
        //    Returns latlon representing Latitude and Longitude.
        
        let Nprime:Double = northing - oN;
        let Mo:Double = calcM(oLat);
        let Mprime:Double = Mo + (Nprime / k);
        let n:Double = (a - b) / (a + b);
        let n2:Double = n * n;
        let n3:Double = n2 * n;
        let n4:Double = n2 * n2;
        let G:Double = a * (1.0 - n) * (1.0 - n2) * (1.0 + (9.0 * n2 / 4.0) + (225.0 * n4 / 64.0)) * (PI / 180.0);
        let sigma:Double = (Mprime * PI) / (180.0 * G);
        
        let latPrimeT1:Double = ((3.0 * n / 2.0) - (27.0 * n3 / 32.0)) * sin(2.0 * sigma);
        let latPrimeT2:Double = ((21.0 * n2 / 16.0) - (55.0 * n4 / 32.0)) * sin(4.0 * sigma);
        let latPrimeT3:Double = (151.0 * n3 / 96.0) * sin(6.0 * sigma);
        let latPrimeT4:Double = (1097.0 * n4 / 512.0) * sin(8.0 * sigma);
        let latPrime:Double = sigma + latPrimeT1 + latPrimeT2 + latPrimeT3 + latPrimeT4;
        
        let sinLatPrime:Double = sin(latPrime);
        let sin2LatPrime:Double = sinLatPrime * sinLatPrime;
        
        let rhoPrime:Double = calcRho(sin2LatPrime);
        let vPrime:Double = calcV(sin2LatPrime);
        let psiPrime:Double = vPrime / rhoPrime;
        let psiPrime2:Double = psiPrime * psiPrime;
        let psiPrime3:Double = psiPrime2 * psiPrime;
        let psiPrime4:Double = psiPrime3 * psiPrime;
        let tPrime:Double = tan(latPrime);
        let tPrime2:Double = tPrime * tPrime;
        let tPrime4:Double = tPrime2 * tPrime2;
        let tPrime6:Double = tPrime4 * tPrime2;
        let Eprime:Double = easting - oE;
        let x:Double = Eprime / (k * vPrime);
        let x2:Double = x * x;
        let x3:Double = x2 * x;
        let x5:Double = x3 * x2;
        let x7:Double = x5 * x2;
        
        // Compute Latitude
        let latFactor:Double = tPrime / (k * rhoPrime);
        let latTerm1:Double = latFactor * ((Eprime * x) / 2.0);
        let latTerm2:Double = latFactor * ((Eprime * x3) / 24.0) * ((-4.0 * psiPrime2) + (9.0 * psiPrime) * (1.0 - tPrime2) + (12.0 * tPrime2));
        let latTerm3:Double = latFactor * ((Eprime * x5) / 720.0) * ((8.0 * psiPrime4) * (11.0 - 24.0 * tPrime2) - (12.0 * psiPrime3) * (21.0 - 71.0 * tPrime2) + (15.0 * psiPrime2) * (15.0 - 98.0 * tPrime2 + 15.0 * tPrime4) + (180.0 * psiPrime) * (5.0 * tPrime2 - 3.0 * tPrime4) + 360.0 * tPrime4);
        let latTerm4:Double = latFactor * ((Eprime * x7) / 40320.0) * (1385.0 - 3633.0 * tPrime2 + 4095.0 * tPrime4 + 1575.0 * tPrime6);
        
        let latRad:Double = latPrime - latTerm1 + latTerm2 - latTerm3 + latTerm4;
        let lat:Double = latRad / (PI / 180.0);
        
        // Compute Longitude
        let secLatPrime:Double = 1.0 / cos(latRad);
        let lonTerm1:Double = x * secLatPrime;
        let lonTerm2:Double = ((x3 * secLatPrime) / 6.0) * (psiPrime + 2.0 * tPrime2);
        let lonTerm3:Double = ((x5 * secLatPrime) / 120.0) * ((-4.0 * psiPrime3) * (1.0 - 6.0 * tPrime2) + psiPrime2 * (9.0 - 68.0 * tPrime2) + 72.0 * psiPrime * tPrime2 + 24.0 * tPrime4);
        let lonTerm4:Double = ((x7 * secLatPrime) / 5040.0) * (61.0 + 662.0 * tPrime2 + 1320.0 * tPrime4 + 720.0 * tPrime6);
        
        let lon:Double = ((oLon * PI / 180.0) + lonTerm1 - lonTerm2 + lonTerm3 - lonTerm4) / (PI / 180.0);
        
        return (lat, lon)
    }
}
