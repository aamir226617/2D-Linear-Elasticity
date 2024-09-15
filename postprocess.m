function postprocess(conM, corM, nodM, sol)
%Material property
E = 3e7 ;
nu=0.3 ;
D = (E / ( 1 - nu^2 ) ) * [1 nu 0; nu 1 0; 0 0 (1-nu) / 2];
% POST-PROCESSING OF THE SOLUTION
% Lets go to each element and prepare a grid in zi and eta
[zi, eta] = meshgrid( -1 : 1 : 1 );
% Row-wise shape fns are arranged
Zi = reshape( zi , [ ],1 );
Eta = reshape( eta , [ ], 1 );
% Obtain values of Ni's at the grid points in master element
N2d = shapeFn2d([ Zi Eta ]);                                    
% Obtain gradient of Ni's at the grid points in master element
dN2d = dshapeFn([ Zi Eta ]);                 
% Number of elements
nElems = size( conM, 1 );
for k = 1 : nElems
    gNods = nodM( k, : );
    xCordNods = corM( gNods, 1);
    yCordNods = corM( gNods, 2);    
    % x and y coordinates corresponding to zi and eta
    xZi = N2d * xCordNods;
    yZi = N2d * yCordNods;    
    % Obtain solution values at those points from (alpha_i* Ni)
    gDofs = conM( k , :);
    alphas = sol( gDofs );
    for p = 1 : length( Zi )
        N = [N2d( p , 1)  0  N2d( p , 2)  0  N2d( p , 3)  0  N2d( p , 4)  0 ; 
                0  N2d( p , 1)  0  N2d( p , 2)  0   N2d( p , 3)  0  N2d( p , 4) ];
        U( : , p) = N * alphas;
    end
    % Surface plot of Temperature
    X = reshape( xZi, size( zi ));
    Y = reshape( yZi, size( zi ));
    Ux = reshape( U( 1 , : ) , size( zi ));
    Uy = reshape( U( 2 , : ) , size( zi ));
    % Plot the solution in this element
%     figure( 3 )
%     hold on
%     surf( X , Y , Ux)
    
    figure( 4 )
    hold on
    surf( X , Y , Uy)
%     % Lets calculate the gradients Strain
    cordXY=[xCordNods yCordNods]';
    for p = 1 : length(Zi)                                                                            
        J = cordXY * dN2d( : , : , p)';                      % Jacobian
        dN = ( J \ eye( 2 ) )' * dN2d( : , : , p);       % derivative in x and y     
        B = [dN( 1 , 1)     0       dN(1 , 2 )      0       dN( 1 , 3 )     0       dN( 1 , 4 )     0; 
                0       dN( 2 , 1 )     0       dN( 2 , 2 )     0       dN( 2 , 3 )     0       dN( 2 , 4 )  ;
                dN( 2 , 1 )  dN( 1 , 1 )  dN( 2 , 2 )  dN( 1, 2 )  dN( 2 , 3 )  dN( 1 , 3 )  ...
                                                                                                    dN( 2 , 4 ) dN( 1 , 4 ) ];
        strn( : , p ) = B * alphas;   
    end
   strs = D * strn;
    strnx = reshape( strn( 1 , : ), size( zi ));
    strny = reshape( strn( 2 , : ), size( zi ));
    strnxy = reshape( strn( 3 , : ), size( zi ));
    
    strsx = reshape( strs( 1 , : ), size( zi ));
    strsy = reshape( strs( 2 , : ), size( zi ));
    strsxy = reshape( strs( 3 , : ), size( zi ));

    figure( 5 )
    hold on
    surf( X , Y , strnx)

%     figure( 6 )
%     hold on
%     surf( X , Y , strny)
%     
%     figure( 7 )
%     hold on 
%     surf( X , Y , strnxy)
    
    figure( 8 )
    hold on
    surf( X , Y , strsx) 
    
%     figure( 9 )
%     hold on
%     surf( X , Y , strsy)
% 
%     figure( 10 )
%     hold on
%     surf( X , Y , strsxy)
end

%         figure( 3 )
%         title( ' Displacement (Ux) ')
%         xlabel( 'X')
%         ylabel( 'Y' )
%         zlabel('Ux ')

         figure( 4 )
          title( ' Displacement (Uy) ')
         xlabel( 'X')
         ylabel( 'Y' )
         zlabel(' Uy ')

        figure( 5 )
        title(' Strain (\epsilon_{xx})')
        xlabel( 'X')
        ylabel( 'Y' )
        zlabel(' \epsilon_{xx}')

%         figure( 6 )
%         title(' Strain (\epsilon_{yy})')
%         xlabel( 'X')
%         ylabel( 'Y' )
%         zlabel(' \epsilon_{yy}')   
% 
%         figure( 7 )
%         title(' Strain (\epsilon_{xy})')
%         xlabel( 'X')
%         ylabel( 'Y' )
%         zlabel(' \epsilon_{xy}')

        figure( 8 )
        title('Stress (\sigma_{xx})')
        xlabel( 'X')
        ylabel( 'Y' )
        zlabel(' \sigma_{xx}')    

%         figure( 9 )
%         title('Stress (\sigma_{yy})')
%         xlabel( 'X')
%         ylabel( 'Y' )
%         zlabel(' \sigma_{yy}')  
% 
%         figure( 10 )
%          title('Stress (\sigma_{xy})')
%         xlabel( 'X')
%         ylabel( 'Y' )
%         zlabel(' \sigma_{xy}')  
    








