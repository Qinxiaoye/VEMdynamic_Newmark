function [point,weight] = gaussInt(nip)
% return the coor for gauss integral point and weight

if nip == 3
    point = [-0.774596669241483,0,0.774596669241483]; % 积分点
    weight = [0.55555555555556,0.888888888888888889,0.55555555555556]; % 权系数
elseif nip == 1
    point = 0;
    weight = 2;
elseif nip == 4
    point = [-0.8611363115940520,-0.3399810435848560,0.3399810435848560,0.8611363115940520];
    weight = [0.3478548451374530,0.6521451548625460,0.6521451548625460,0.3478548451374530];
elseif nip == 5
    point = [-0.9061798459386640;-0.5384693101056830;0.0000000000000000;0.5384693101056830;0.9061798459386640];
    weight = [0.2369268850561890;0.4786286704993660;0.5688888888888880;0.4786286704993660;0.2369268850561890];
elseif nip == 6
    point = [-0.9324695142031520;-0.6612093864662640;-0.2386191860831960;0.2386191860831960;0.6612093864662640;0.9324695142031520];
    weight = [0.1713244923791700;0.3607615730481380;0.4679139345726910;0.4679139345726910;0.3607615730481380;0.1713244923791700];
elseif nip == 7
    point = [-0.9491079123427580;-0.7415311855993940;-0.4058451513773970;0.0000000000000000;0.4058451513773970;0.7415311855993940;0.9491079123427580];
    weight = [0.1294849661688690;0.2797053914892760;0.3818300505051180;0.4179591836734690;0.3818300505051180;0.2797053914892760;0.1294849661688690];
elseif nip == 8
    point = [-0.9602898564975360;-0.7966664774136260;-0.5255324099163290;-0.1834346424956490;0.1834346424956490;0.5255324099163290;0.7966664774136260;0.9602898564975360];
    weight = [0.1012285362903760;0.2223810344533740;0.3137066458778870;0.3626837833783620;0.3626837833783620;0.3137066458778870;0.2223810344533740;0.1012285362903760];
elseif nip == 9
    point = [-0.9681602395076260;-0.8360311073266350;-0.6133714327005900;-0.3242534234038080;0.0000000000000000;0.3242534234038080;0.6133714327005900;0.8360311073266350;0.9681602395076260];
    weight = [0.0812743883615744;0.1806481606948570;0.2606106964029350;0.3123470770400020;0.3302393550012590;0.3123470770400020;0.2606106964029350;0.1806481606948570;0.0812743883615744];


end