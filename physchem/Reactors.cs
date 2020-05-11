using System;

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
    abstract class CourseReactor : IReactor
    {
        public event EventHandler ProgressCallback;

        public double[] AdsorbateConcentration { get; set; }
        public double[] OccupiedSiteFraction { get; set; }
        public double[] TimeStamps { get; set; }

        public double CurrentTime { get; protected set; }
        public double TimeStep = 0.001;

        public virtual void CalculateRegion(double duration)
        {
            int length = Initialize(duration) - 1;
            SetInitialConditions();
            for (int i = 0; i < length; i++)
            {
                PerformStep(i, TimeStep);
            }
        }

        protected virtual int Initialize(double duration)
        {
            int length = (int)Math.Floor(duration / TimeStep);
            AdsorbateConcentration = new double[length];
            OccupiedSiteFraction = new double[length];
            TimeStamps = new double[length];
            TimeStamps[0] = 0;  
            CurrentTime = 0;
            return length;
        }
        protected abstract void PerformStep(int i, double timeStep);
        protected abstract void SetInitialConditions();
    }

    class ReactorForDesorption : CourseReactor
    {
        public double DebyeLength;
        public double BulkElectronConcentration;
        public double AdsorptionConstant;
        public double DesorptionConstant;
        public double VolumeFlowRate;
        public double SurfaceArea;
        public double ReactorVolume;

        public double InitialOccupiedSitesFraction;
        public double InitialConcentrationOfSites;

        protected override void SetInitialConditions()
        {
            AdsorbateConcentration[0] = 0;
            OccupiedSiteFraction[0] = InitialOccupiedSitesFraction;
        }

        protected override void PerformStep(int i, double timeStep)
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

    class ReactorForPhotoDesorption : ReactorForDesorption
    {
        //public double OxygenVacancyConcentration;
        public double GenerationProduct;

        protected override void PerformStep(int i, double timeStep)
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

    class ReactorForAdsorption : CourseReactor
    {
        public double DebyeLength;
        public double BulkElectronConcentration;
        public double AdsorptionConstant;
        public double VolumeFlowRate;
        public double SurfaceArea;
        public double ReactorVolume;

        public double InitialOccupiedSitesFraction;
        public double InitialConcentrationOfSites;
        public double InitialAdsorbateConcentration;

        protected override void SetInitialConditions()
        {
            AdsorbateConcentration[0] = InitialAdsorbateConcentration;
            OccupiedSiteFraction[0] = InitialOccupiedSitesFraction;
        }

        protected override void PerformStep(int i, double timeStep)
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

    class ReactorForPhotoAdsorption : ReactorForAdsorption
    {
        public double GenerationProduct;

        protected override void PerformStep(int i, double timeStep)
        {
            //Compute A'(t)
            double aPrime = AdsorptionConstant * AdsorbateConcentration[i] * (1 - OccupiedSiteFraction[i]) *
                (BulkElectronConcentration + GenerationProduct) *
                Math.Exp(-InitialConcentrationOfSites * OccupiedSiteFraction[i] / (DebyeLength * BulkElectronConcentration));
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
}
