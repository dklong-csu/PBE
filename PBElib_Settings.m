function settings = PBElib_Settings(vecSize, options)
%NanoSimSettings settings to control particle simulation
    arguments
        vecSize {mustBeInteger} = 1
        options.A {mustBeReal} = zeros(vecSize,1)
        options.K {mustBeReal} = 0
        
        options.pstart {mustBeInteger,mustBePositive} = 1
        options.pend {mustBeInteger,mustBePositive} = 1
        options.pBins {mustBeVector,mustBeReal} = [1,1]
        options.gidx {mustBeScalarOrEmpty,mustBeInteger,mustBePositive} = []
        options.gKernel {mustBeReal} = []
        options.gRxnIdx {mustBeInteger} = []
        options.gRxnCoeff {mustBeInteger} = []
        options.sizes {mustBeFloat,mustBeNonnegative} = [];
        options.dsizes {mustBeFloat,mustBeNonnegative} = [];

        options.cutoff {mustBeScalarOrEmpty, mustBeInteger} = []
        options.aKernel {mustBeReal} = []
        options.aMapping cell = {}

        options.useEMoM {mustBeNumericOrLogical, mustBeScalarOrEmpty} = false
        options.emomDelx {mustBeFloat,mustBeFinite,mustBePositive} = 1.0
        options.emomInflowRate {mustBeFloat, mustBeFinite,mustBeNonnegative} = 1.0
        options.emomGrowthRate {mustBeFloat, mustBeFinite, mustBePositive} = 1.0
        options.emomBigParticleDiam {mustBeFloat, mustBeFinite, mustBePositive} = 1.0
    end

    settings.vecSize = vecSize;
    settings.A = options.A;
    settings.K = options.K;

    settings.pstart = options.pstart;
    settings.pend = options.pend;
    settings.pBins = options.pBins;
    settings.gidx = options.gidx;
    settings.gKernel = options.gKernel;
    settings.gRxnIdx = options.gRxnIdx;
    settings.gRxnCoeff = options.gRxnCoeff;
    settings.sizes = options.sizes;
    settings.dsizes = options.dsizes;

    settings.cutoff = options.cutoff;
    settings.aKernel = options.aKernel;
    settings.aMapping = options.aMapping;

    settings.useEMoM = options.useEMoM;
    settings.emomDelx = options.emomDelx;
    settings.emomInflowRate = options.emomInflowRate;
    settings.emomGrowthRate = options.emomGrowthRate;
    settings.emomBigParticleDiam = options.emomBigParticleDiam;
end
