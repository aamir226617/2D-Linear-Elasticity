function LocalNodsOnEdge=ElementLocalNumberingScheme(degP)
% Local nodes on edges
switch degP
    case 1
        LocalNodsOnEdge=[1  2;  
                                             2   3;  
                                            3   4;  
                                            4   1];
    case 2
        LocalNodsOnEdge=[1  2   3;
                                            3   4   5;
                                            5   6   7;
                                            7   8   1];
end