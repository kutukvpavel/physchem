#define DEBUG

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace physchem
{
    interface IReactor
    {
        event EventHandler ProgressCallback;
        void CalculateRegion(double duration);
        double[] AdsorbateConcentration { get; }
        double[] OccupiedSiteFraction { get; }
        double[] TimeStamps { get; }
    }

    class ReactorForPhotoAdsorption : ReactorForAdsorption, IReactor
    {
        public double GenerationProduct;

        public void CalculateRegion(double duration)
        {
            int length = (int)Math.Floor(duration / TimeStep);
            AdsorbateConcentration = new double[length];
            OccupiedSiteFraction = new double[length];
            TimeStamps = new double[length];
            TimeStamps[0] = 0;
            AdsorbateConcentration[0] = InitialAdsorbateConcentration;
            OccupiedSiteFraction[0] = InitialOccupiedSitesFraction;
            CurrentTime = 0;
            length--;
            for (int i = 0; i < length; i++)
            {
                PerformStep(i, TimeStep);
            }
        }

        public void PerformStep(int i, double timeStep)
        {
            //Compute A'(t)
            double aPrime = AdsorptionConstant * AdsorbateConcentration[i] * (1 - OccupiedSiteFraction[i]) *
                (BulkElectronConcentration + GenerationProduct) * 
                Math.Exp(-InitialOccupiedSitesFraction * OccupiedSiteFraction[i] / (DebyeLength * BulkElectronConcentration));
            //Compute C'(t)
            double cPrime = -InitialConcentrationOfSites * SurfaceArea * aPrime +
                VolumeFlowRate * (InitialAdsorbateConcentration - AdsorbateConcentration[i]) / ReactorVolume;
            //Compute linear change
            AdsorbateConcentration[i + 1] = AdsorbateConcentration[i] + cPrime * timeStep;
            OccupiedSiteFraction[i + 1] = OccupiedSiteFraction[i] + aPrime * timeStep;
            //Compute time
            CurrentTime += timeStep;
            TimeStamps[i + 1] = CurrentTime;
        }
    }

    class ReactorForDesorption : IReactor
    {
        public event EventHandler ProgressCallback;

        public double[] AdsorbateConcentration { get; set; }
        public double[] OccupiedSiteFraction { get; set; }
        public double[] TimeStamps { get; set; }

        public double DebyeLength;
        public double BulkElectronConcentration;
        public double AdsorptionConstant;
        public double DesorptionConstant;
        public double VolumeFlowRate;
        public double SurfaceArea;
        public double ReactorVolume;

        public double InitialOccupiedSitesFraction;
        public double InitialConcentrationOfSites; 

        public double CurrentTime;
        public double TimeStep = 0.001;

        public void CalculateRegion(double duration)
        {
            int length = (int)Math.Floor(duration / TimeStep);
            AdsorbateConcentration = new double[length];
            OccupiedSiteFraction = new double[length];
            TimeStamps = new double[length];
            TimeStamps[0] = 0;
            AdsorbateConcentration[0] = 0;
            OccupiedSiteFraction[0] = InitialOccupiedSitesFraction;
            CurrentTime = 0;
            length--;
            for (int i = 0; i < length; i++)
            {
                PerformStep(i, TimeStep);
            }
        }

        public void PerformStep(int i, double timeStep)
        {
            //Compute A'(t)
            double aPrime = AdsorptionConstant * AdsorbateConcentration[i] * (1 - OccupiedSiteFraction[i]) *
                (BulkElectronConcentration - OccupiedSiteFraction[i] * InitialConcentrationOfSites / DebyeLength) -
                DesorptionConstant * OccupiedSiteFraction[i] *
                (BulkElectronConcentration + OccupiedSiteFraction[i] * InitialConcentrationOfSites / DebyeLength);
            //Compute C'(t)
            double cPrime = -InitialConcentrationOfSites * SurfaceArea * aPrime -
                VolumeFlowRate * AdsorbateConcentration[i] / ReactorVolume;
            //Compute linear change
            AdsorbateConcentration[i + 1] = AdsorbateConcentration[i] + cPrime * timeStep;
            OccupiedSiteFraction[i + 1] = OccupiedSiteFraction[i] + aPrime * timeStep;
            //Compute time
            CurrentTime += timeStep;
            TimeStamps[i + 1] = CurrentTime;
        }
    }

    class ReactorForAdsorption : IReactor
    {
        public event EventHandler ProgressCallback;

        public double[] AdsorbateConcentration { get; set; }
        public double[] OccupiedSiteFraction { get; set; }
        public double[] TimeStamps { get; set; }

        public double DebyeLength;
        public double BulkElectronConcentration;
        public double AdsorptionConstant;
        public double VolumeFlowRate;
        public double SurfaceArea;
        public double ReactorVolume;

        public double InitialOccupiedSitesFraction;
        public double InitialConcentrationOfSites;
        public double InitialAdsorbateConcentration;

        public double CurrentTime;
        public double TimeStep = 0.001;

        public void CalculateRegion(double duration)
        {
            int length = (int)Math.Floor(duration / TimeStep);
            AdsorbateConcentration = new double[length];
            OccupiedSiteFraction = new double[length];
            TimeStamps = new double[length];
            TimeStamps[0] = 0;
            AdsorbateConcentration[0] = InitialAdsorbateConcentration;
            OccupiedSiteFraction[0] = InitialOccupiedSitesFraction;
            CurrentTime = 0;
            length--;
            for (int i = 0; i < length; i++)
            {
                PerformStep(i, TimeStep);
            }
        }

        public void PerformStep(int i, double timeStep)
        {
            //Compute A'(t)
            double aPrime = AdsorptionConstant * AdsorbateConcentration[i] * (1 - OccupiedSiteFraction[i]) *
                (BulkElectronConcentration - OccupiedSiteFraction[i] * InitialConcentrationOfSites / DebyeLength);
            //Compute C'(t)
            double cPrime = -InitialConcentrationOfSites * SurfaceArea * aPrime + 
                VolumeFlowRate * (InitialAdsorbateConcentration - AdsorbateConcentration[i]) / ReactorVolume;
            //Compute linear change
            AdsorbateConcentration[i + 1] = AdsorbateConcentration[i] + cPrime * timeStep;
            OccupiedSiteFraction[i + 1] = OccupiedSiteFraction[i] + aPrime * timeStep;
            //Compute time
            CurrentTime += timeStep;
            TimeStamps[i + 1] = CurrentTime;
        }
    }

    class ReactorForPhotoDesorption : IReactor
    {
        public event EventHandler ProgressCallback;

        public double[] AdsorbateConcentration { get; set; }
        public double[] OccupiedSiteFraction { get; set; }
        public double[] TimeStamps { get; set; }

        public double DebyeLength;
        public double BulkElectronConcentration;
        //public double OxygenVacancyConcentration;
        public double GenerationProduct;
        public double AdsorptionConstant;
        public double DesorptionConstant;
        public double VolumeFlowRate;
        public double SurfaceArea;
        public double ReactorVolume;

        public double InitialOccupiedSitesFraction;
        public double InitialConcentrationOfSites;

        public double CurrentTime;
        public double TimeStep = 0.001;

        public void CalculateRegion(double duration)
        {
            int length = (int)Math.Floor(duration / TimeStep);
            AdsorbateConcentration = new double[length];
            OccupiedSiteFraction = new double[length];
            TimeStamps = new double[length];
            TimeStamps[0] = 0;
            AdsorbateConcentration[0] = 0;
            OccupiedSiteFraction[0] = InitialOccupiedSitesFraction;
            CurrentTime = 0;
            length--;
            for (int i = 0; i < length; i++)
            {
                PerformStep(i, TimeStep);
            }
        }

        public void PerformStep(int i, double timeStep)
        {
            double exponent =
                Math.Exp(-InitialConcentrationOfSites * OccupiedSiteFraction[i] /
                (DebyeLength * BulkElectronConcentration));
            //Compute A'(t)
            double aPrime = AdsorptionConstant * AdsorbateConcentration[i] * (1 - OccupiedSiteFraction[i]) 
                * BulkElectronConcentration * (exponent) -
                DesorptionConstant * OccupiedSiteFraction[i] * GenerationProduct * (1 / exponent);
            //Compute C'(t)
            double cPrime = -InitialConcentrationOfSites * SurfaceArea * aPrime -
                VolumeFlowRate * AdsorbateConcentration[i] / ReactorVolume; 
            /*double exponent = InitialConcentrationOfSites * OccupiedSiteFraction[i] /
                (DebyeLength * BulkElectronConcentration);
            //Compute A'(t)
            double aPrime = AdsorptionConstant * AdsorbateConcentration[i] * (1 - OccupiedSiteFraction[i])
                * BulkElectronConcentration * (1 - exponent) -
                DesorptionConstant * OccupiedSiteFraction[i] * GenerationProduct * (1 + exponent);
            //Compute C'(t)
            double cPrime = -InitialConcentrationOfSites * SurfaceArea * aPrime -
                VolumeFlowRate * AdsorbateConcentration[i] / ReactorVolume;  */
            //Compute linear change
            AdsorbateConcentration[i + 1] = AdsorbateConcentration[i] + cPrime * timeStep;
            OccupiedSiteFraction[i + 1] = OccupiedSiteFraction[i] + aPrime * timeStep;
            //Compute time
            CurrentTime += timeStep;
            TimeStamps[i + 1] = CurrentTime;
        }
    }

    class Program
    {
        public const double AvogadroNumber = 6.02E23;

        static ReactorForAdsorption r = new ReactorForAdsorption()
        {
            InitialAdsorbateConcentration = 0.0044,       //mol/m^3 O2 100ppm
            InitialConcentrationOfSites = 1E18 / AvogadroNumber, //mol/m^2
            SurfaceArea = 300, //m^2/sample
            InitialOccupiedSitesFraction = 0,
            DebyeLength = 20E-9,                   //m
            BulkElectronConcentration = 1E18 * 1E6 / AvogadroNumber,   //mol/m^-3
            AdsorptionConstant = 1,          //s^-1
            VolumeFlowRate = 10E-6 / 60,          //m^3/s
            TimeStep = 1E-4,
            ReactorVolume = 32E-6 //m^3
        };

        static ReactorForDesorption r2 = new ReactorForDesorption()
        {
            InitialConcentrationOfSites = 1E18 / AvogadroNumber, //mol/m^2
            SurfaceArea = 300, //m^2/sample
            InitialOccupiedSitesFraction = 0.8,
            DebyeLength = 20E-9,                   //m
            BulkElectronConcentration = 1E18 * 1E6 / AvogadroNumber,   //mol/m^-3
            AdsorptionConstant = 1,          //s^-1
            DesorptionConstant = 0.01,
            VolumeFlowRate = 10E-6 / 60,          //m^3/s
            TimeStep = 1E-4,
            ReactorVolume = 32E-6 //m^3
        };

        static ReactorForPhotoDesorption r3 = new ReactorForPhotoDesorption()
        {
            InitialConcentrationOfSites = 7E10 / AvogadroNumber, //mol/m^2
            SurfaceArea = 300, //m^2/sample
            InitialOccupiedSitesFraction = 0.8,
            DebyeLength = 2E-6,                   //m
            BulkElectronConcentration = 1E13 * 1E6 / AvogadroNumber,   //mol/m^-3
            AdsorptionConstant = 1E6,          //s^-1
            DesorptionConstant = 1E6,
            VolumeFlowRate = 10E-6 / 60,          //m^3/s
            TimeStep = 1E-4,
            ReactorVolume = 32E-6, //m^3
            GenerationProduct = 1E11 * 1E6 / AvogadroNumber
        };

        static void Main(string[] args)
        {
            Console.WriteLine("Computing...");

            Compute(r3, 100, 100);

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
