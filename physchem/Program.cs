using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace physchem
{
    class Program
    {
        public const double AvogadroNumber = 6.02E23;

        static ReactorForPhotoDesorption r3 = new ReactorForPhotoDesorption()
        {
            InitialConcentrationOfSites = 7E10 / AvogadroNumber, //mol/m^2
            SurfaceArea = 300, //m^2/sample
            InitialOccupiedSitesFraction = 0.8,
            DebyeLength = 2E-6,                   //m
            BulkElectronConcentration = 1E13 * 1E6 / AvogadroNumber,   //mol/m^-3
            AdsorptionConstant = 1E6,          
            DesorptionConstant = 1E6,
            VolumeFlowRate = 10E-6 / 60,          //m^3/s
            TimeStep = 1E-4,
            ReactorVolume = 32E-6, //m^3
            GenerationProduct = 1E11 * 1E6 / AvogadroNumber
        };
        static ReactorForPhotoAdsorption r4 = new ReactorForPhotoAdsorption()
        {
            InitialConcentrationOfSites = 7E10 / AvogadroNumber, //mol/m^2
            SurfaceArea = 300, //m^2/sample
            InitialOccupiedSitesFraction = 0,
            DebyeLength = 2E-6,                   //m
            BulkElectronConcentration = 1E13 * 1E6 / AvogadroNumber,   //mol/m^-3
            AdsorptionConstant = 1E9,          
            VolumeFlowRate = 10E-6 / 60,          //m^3/s
            TimeStep = 1E-4,
            ReactorVolume = 32E-6, //m^3
            GenerationProduct = 1E11 * 1E6 / AvogadroNumber,
            InitialAdsorbateConcentration = 1E-6 //mol/m^3
        };

        static void Main(string[] args)
        {
            Console.WriteLine("Computing...");

            Compute(r4, 300, 100);

            Console.WriteLine("Done.");
        }

        static void Compute(IReactor r, int seconds, int stepDivisor)
        {
            r.CalculateRegion(seconds);
            StringBuilder results = new StringBuilder();
            int skipped = 0;
            for (int i = 0; i < r.TimeStamps.Length; i++)
            {
                if (i % stepDivisor == 0) skipped++;
                if (skipped > 0 && skipped < stepDivisor)
                {
                    skipped++;
                }
                else
                {
                    skipped = 0;
                    results.AppendFormat("{0:F4};{1:F7};{2:F6}" + Environment.NewLine,
                        r.TimeStamps[i], r.AdsorbateConcentration[i] * 1E6, r.OccupiedSiteFraction[i]);
                }
            }
            File.WriteAllText(@"E:\course.csv", results.ToString());
        }
        
    }
}
