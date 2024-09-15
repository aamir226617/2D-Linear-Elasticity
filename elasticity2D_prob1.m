clear all
close all
clc
% FE CODE FOR SOLVING POISSON'S EQUATION IN 2D
% IT USES QUADRILATERAL ELEMENT

%Material property
E = 210e9 ;
nu=0.3 ;
D = (E / ( 1 - nu^2 ) ) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
%Body force term
b0 = @(x,y) [0.0; 0.0];
% Define values of Neumann boundary     (Only force along X-axis on right side, others are homogeneous NB)
theta=atan(4);
Tn=100;
% bcNv = {@(x,y) [0; 0], @(x,y) [-Tn * sin(theta); -Tn * cos(theta)], @(x,y) [0; 0]};
bcNv = {@(x,y) [0; 0], @(x,y) [360e6; 0], @(x,y) [0; 0]};
%  Value of u in the Dirichlet boundary
bcDv = 0; 

% No. of nodes
degP = 1;
nNod_elem = 4;                            % no. of nodes per element
dof_Nod = 2;                                     % no. of dof per node
ndof_elem = nNod_elem * dof_Nod;     % no. of dofs per element
% Local node numbering scheme
localNodsOnEdge = ElementLocalNumberingScheme( degP );                                     
 
% Call refine and get the mesh data
nRefineStep = 4;                                            % refine it two times
refine( nRefineStep , degP, dof_Nod);

% Gauss Point in 2D
gp2d = [ -1/sqrt(3.0)      -1/sqrt(3.0) ;
                1/sqrt(3.0)       -1/sqrt(3.0) ;
                1/sqrt(3.0)        1/sqrt(3.0) ;
               -1/sqrt(3.0)        1/sqrt(3.0) ];
 % Corresponding weight
 w2d = [   1.0;
                 1.0;
                 1.0;
                 1.0 ]; 
 % Gauss Point in 1D
 gp1d = [-1/sqrt(3.0)        1/sqrt(3.0)]';
 % Weights in 1D
 w1d = [1 1]';      
 
 %Evaluate shape fns at Gauss pts
N2d = shapeFn2d( gp2d );   
%Evaluate derivative of the Shape fns
dN2d = dshapeFn( gp2d ) ;
%Evaluate 1D shape fns at Gauss points
N1d = shapeFn1d( gp1d );

for refine_step=1 : nRefineStep
    refine_step
    fname = sprintf( 'meshInfo%d.mat' , refine_step );
    load( fname );                                                           % loads the mesh information conM, corM, bcTyp
    nElems = size( conM , 1 );
    %Total nodes and DOF
    nNods = size( corM , 1 );
    nDofs = nNods * dof_Nod;

    % Global stiffness matrix and force vector
    kg = zeros( nDofs );
    fg = zeros( nDofs , 1);  
  
for k =1 : nElems
     ke = zeros( ndof_elem );
     fe = zeros(ndof_elem,1);
     cornerNods = localNodsOnEdge( : , 1 );
     % X and Y coordinates of the corner nodes (LINEAR MAP)
     cordXY = corM( nodM( k, cornerNods), : );
     % FOR EACH GAUSS POINT
     for p = 1 : size( gp2d , 1 )          
         % Physical Y-Cords for GaussPoints
         xyGp = N2d( p, : ) * cordXY;
         % Calculating Lambda0 & source term g0 at (xp,yp)
         b = b0( xyGp(1) , xyGp(2));         
         % Calculating Jacobian & its determinant at the Gauss Point
         J = cordXY' * dN2d( : , : , p)';
         detJ = det( J );
         % Derivatives of shape fns in physical element
         dN = (J \ eye( 2 ))' * dN2d( : , : , p);
         % Now construct B matrix
         B = [dN( 1 , 1)      0                   dN(1 , 2 )            0                dN( 1 , 3 )       0                     dN( 1 , 4 )     0                 ;
                0                    dN( 2 , 1 )     0                       dN( 2 , 2 )    0                     dN( 2 , 3 )      0                    dN( 2 , 4 )  ;
                dN( 2 , 1 )     dN( 1 , 1 )     dN( 2 , 2 )         dN( 1, 2 )    dN( 2 , 3 )       dN( 1 , 3 )      dN( 2 , 4 )     dN( 1 , 4 )  ];
         % Construct N matrix
         N = [N2d( p , 1)   0                   N2d( p , 2)     0                   N2d( p , 3)     0                   N2d( p , 4)     0                  ;
             0                  N2d( p , 1)   0                     N2d( p , 2)     0                   N2d( p , 3)     0                   N2d( p , 4) ];
         
         % CALCULATING ELEMENT STIFFNESS AND FORCE ENTRIES
         ke = ke + B' * D * B * detJ * w2d( p );
         fe = fe + N' * b * detJ * w2d( p );
     end
     % CALCULATING THE CONTRIBUTION OF NEUMANN BC (IF ANY)
     % finding all edges on Neumann boundary having bcTyp positive
     nLocalEdge = find( bcTyp( k , : ) > 0);                                                                                % nLocal_edge ---->> local edge numbers  
     % if atleast one edge coincides with Neumann boundary
     if (~isempty( nLocalEdge ))         
            % loop over all such edges
             for j = 1 : length( nLocalEdge )         
                 % global nodes on the j^th edge                  
                globalNodsOnEdge = nodM( k, localNodsOnEdge( nLocalEdge( j ), : ));         
                % coordinates of the end nodes of the j^th edge
                XY = corM(globalNodsOnEdge( [1 end] ), :);
                % length of the j^th edge
                h = sqrt( (diff( XY( : , 1) ))^2 + (diff( XY( : , 2)))^2 ) /2;
                % loop over 1D Gauss points
                for p = 1 : length( gp1d )
                         xyGp = N1d( p, :) * XY;                % cordinates on physical boundary
                         % evaluate traction T(xp,yp)
                         T = bcNv{ bcTyp( k , nLocalEdge( j ) )} (xyGp (1) , xyGp(2));   
                         % corresponding 1D shape fns at p^th Gauss points
                         N1 = zeros(dof_Nod, ndof_elem);
                         local_nods = localNodsOnEdge( nLocalEdge( j ), : );
                         N1(1, 2*local_nods-1) = N1d(p, : );
                         N1(2, 2*local_nods) = N1d(p, : );
                         % calculate NB contribution using Gauss quadrature
                        fe = fe + N1' * T * h * w1d( p );
                end
            end
     end
% ke
% fe
    % ASSEMBLY
     for i = 1 : ndof_elem
            ig = conM( k , i );
            for j = 1 : ndof_elem
                    jg = conM( k , j );
                    kg( ig , jg ) = kg( ig , jg ) + ke( i , j);
             end
             fg( ig ) = fg( ig ) + fe( i );
     end     
end

% % APPLY DIRICHLET BC
% [el, ed] = find( bcTyp < 0 );        % element and edge on Dirichlet BC
% nodDb = unique( nodM( el , localNodsOnEdge( ed , : ) ) );
% dofDb = reshape( [(2 * nodDb -1)     (2 * nodDb)] , [  ] , 1);
dofDb = [1 2 4 7];
for i = 1 : length( dofDb )
         j = dofDb( i );
         fg = fg - bcDv * kg( : , j );
         kg( j , : ) = 0;
         kg( : ,  j ) = 0;
         kg( j , j ) = 1;
         fg( j ) = bcDv;
end
kg
 %SOLVE THE SYSTEM OF EQUATION
 fg
alpha = kg \ fg;
 alpha(3)
 U(refine_step) = (alpha' * kg * alpha)/2;
 nalpha(refine_step) = length( alpha );
end
% 
% % LETS PLOT THE SOLUTION
% 
figure(2)
plot(nalpha , U , 'o-')
set(findall(gcf,'-property','FontSize'),'FontSize', 18)
xlabel( 'No. of unknowns' )
ylabel( 'U' )
% 
%  
postprocess(conM, corM, nodM, alpha)
% 
% postprocess_gausspt(conM, corM, nodM, alpha)

% VTK_Plane_Elasticity(nodM, corM, conM, alpha)