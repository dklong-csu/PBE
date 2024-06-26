function [a,b] = adhesionModel(r,a0,b0,type)
    %   Numerically calculated for small sizes
    % coverage_fit = [1	0.900000000000000	0.800000000000000	0.700000000000000	0.600000000000000	0.500000000000000	0.400000000000000	0.300000000000000	0.200000000000000	0.100000000000000	0.0100000000000000	0.00100000000000000	0.000100000000000000	1.00000000000000e-05];
    % a_fit = [0.00622808784368107	0.0195847629616278	0.0573948199339693	0.117063079934492	0.199401641198954	0.300357000273945	0.412953209698401	0.530401195093301	0.652711702456829	0.796311613368275	0.975552785265214	0.997527259429598	0.999752508252221	0.999975248704746];
    % b_fit = [-1.30195733601944	-1.01553857665512	-0.997242430261401	-0.892390683273845	-0.778970841222298	-0.654596232089080	-0.517788449933337	-0.373060138860471	-0.231545957417329	-0.106208064668749	-0.0100339026727572	-0.00100058010878498	-0.000100032243084992	-1.00029684091124e-05];
    % 
    % ainterp = griddedInterpolant(flip(coverage_fit),flip(a_fit),'makima','linear');
    % a = min(1,max(0,ainterp(coverage)));
    % 
    % binterp = griddedInterpolant(flip(coverage_fit), flip(b_fit),'makima','linear');
    % b = min(0,binterp(coverage));
    % a = 0.6;

    switch type
        case 1
            smoothstep = r;
        case 2
            smoothstep = 3*r^2 - 2*r^3;
        case 3
            smoothstep = 6*r^5 - 15*r^4 + 10*r^3;
        case 4
            smoothstep = r.^2;
    end
    smoothstep = max(min(1,smoothstep),0);


    a= 1 - a0*smoothstep;
    % a = a0*(1-sqrt(coverage));
    % b = b0*sqrt(coverage);
    b = b0*smoothstep;

end