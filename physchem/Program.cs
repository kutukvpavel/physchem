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
        static ReactorForNumericAdsorption rNum = new ReactorForNumericAdsorption()
        {
            InitialConcentrationOfSites = 8.5E-6, //mol/m^2
            SurfaceArea = 6.4, //m^2/sample
            InitialOccupiedSitesFraction = 0,
            VolumeFlowRate = 15E-6 / 60,          //m^3/s
            TimeStep = 1E-1,
            ReactorVolume = 32E-6, //m^3
            InitialAdsorbateConcentration = 0.0044 //mol/m^3
        };

        static void Main(string[] args)
        {
            Console.WriteLine("Computing...");

            //Compute(r4, 300, 100);
            ComputeFromNumericAdsorbateConc(rNum, 9620, 47500 - 9620, 100, 3.86E-4, 2.46E+1);

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
            File.WriteAllText(Path.Combine(Environment.CurrentDirectory, @"course.csv"), results.ToString());
        }
        static void ComputeFromNumericAdsorbateConc(ReactorForNumericAdsorption r, int startSeconds, int seconds, int stepDivisor, double calK, double calB)
        {
            //Load MS data
            List<double> time = new List<double>(startSeconds + seconds);
            List<double> conc = new List<double>(startSeconds + seconds);
            var lines = File.ReadAllLines(@"C:\Users\Photo\Desktop\EXPORTS\TiO2-photoads-O2\319_TiO2_02_05_24.csv").Skip(1);
            foreach (var item in lines)
            {
                string[] cells = item.Split(';');
                time.Add(double.Parse(cells[0]));
                conc.Add((double.Parse(cells[1]) - calB) * calK);
            }

            //Perform linear interpolation
            //List<double> interpolatedTime = new List<double>(seconds);
            List<double> interpolatedConc = new List<double>(seconds);
            int stopSeconds = startSeconds + seconds;
            for (int i = 0; i < time.Count - 1; i++)
            {
                if (time[i] < startSeconds) continue;
                if (time[i] > stopSeconds) break;
                double interpolationSteps = (time[i + 1] - time[i]) / r.TimeStep;
                for (int j = 0; j < interpolationSteps; j++)
                {
                    //double interpolatedTime = time[i] + internalTimeStep * j;
                    //interpolatedTime.Add();
                    interpolatedConc.Add(conc[i] + j * (conc[i + 1] - conc[i]) / interpolationSteps);
                }
            }

            //Correct baseline drift (2-point approximation for a photoads. experiment)
            double baselineStep = (interpolatedConc[0] - interpolatedConc.Last()) / interpolatedConc.Count;
            for (int i = 0; i < interpolatedConc.Count; i++)
            {
                interpolatedConc[i] += baselineStep * i;
            }

            //Load data into reactor
            r.AdsorbateConcentration = interpolatedConc.ToArray();

            //Compute
            Compute(r, seconds, stepDivisor);
        }
        
    }
}
