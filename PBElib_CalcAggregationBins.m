function [idx,weight] = PBElib_CalcAggregationBins(i,j,pBins,idxBin1) %#codegen
%   FIXME This is slower than it should be...
    %----------------------------------------------------------------------
    %   Get the bin limits
    %----------------------------------------------------------------------
    binI = pBins(i:i+1);
    binJ = pBins(j:j+1);
    %   How many are in each bin?
    lenI = diff(binI);
    lenJ = diff(binJ);
    %   What is the smallest bin size?
    L = min(lenI,lenJ);

    %   What is the minimum created size?
    minsum = binI(1)+binJ(1);
    %   What is the maximum created size?
    maxsum = binI(2)-1 + binJ(2)-1;

    %----------------------------------------------------------------------
    %   Number of ways to create a size follows the pattern:
    %   minsum      *
    %   minsum+1    **
    %   minsum+2    ***
    %   ...         
    %   minsum+L-2  ***...*     <- L-1 *'s
    %   minsum+L-1  ***...**    <- L *'s
    %   minsum+L    ***...**
    %   ...
    %   maxsum-L    ***...**
    %   maxsum-L+1  ***...**
    %   maxsum-L+2  ***...*     <- L-1 *'s
    %   ...
    %   maxsum-2    ***
    %   maxsum-1    **
    %   maxsum      *
    up = [minsum, minsum+L-2];
    down = [maxsum-L+2, maxsum];
    mid = [minsum+L-1,maxsum-L+1];

    newBins = zeros(length(pBins)-1,1);

    %   Can start looking at the larger bin
    firstbin = find(minsum<pBins);
    firstbin = firstbin(1)-1;
    lastbin = find(maxsum<pBins);
    lastbin = lastbin(1)-1;
    for b=firstbin:lastbin
        %   Current bin being checked
        bin = pBins(b:b+1);

        %   What combinations are on the up-slope in this bin?
        binup = [max(bin(1), up(1))
                 min(bin(2)-1, up(2))];
        if binup(2) >= binup(1)
            %   In this case part of the up-slope is in this bin
            %   Count number of ways to get to sum
            M = binup(1) - minsum;
            N = binup(2) - minsum + 1;
            %   Sum of first N integers minus sum of first M integers
            Msum = M * (M+1) / 2;
            Nsum = N * (N+1) / 2;
            newBins(b) = newBins(b) + Nsum - Msum;
        end

        %   What combinations are on the level set in the middle?
        binmid = [max(bin(1), mid(1))
                  min(bin(2)-1, mid(2))];
        if binmid(2) >= binmid(1)
            %   Then there are combination on the level set in this bin
            N = binmid(2) - binmid(1) + 1; % How many
            newBins(b) = newBins(b) + N * L; % Each has L ways to form
        end

        %   What combinations are on the down-slope in this bin?
        bindown = [max(bin(1), down(1))
                   min(bin(2)-1, down(2))];
        if bindown(2) >= bindown(1)
            %   Then there are combinations on the level set in this bin
            %   Count number of ways to get to sum
            M = -(bindown(2) - maxsum);
            N = -(bindown(1) - maxsum - 1);
            Msum = M * (M+1) / 2;
            Nsum = N * (N+1) / 2;
            newBins(b) = newBins(b) + Nsum - Msum;
        end
    end

    %   newBins contains # ways to create a particle in each bin
    %   Total combinations = lenI * lenJ
    %   Percentage = newBins / (lenI * lenJ)
    %   We only need the non-zero ones
    idx = find(newBins > 0) + idxBin1 -1;
    weight = newBins(idx);
    weight = weight / lenI / lenJ;

end