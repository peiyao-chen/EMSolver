function test_chainy

% Admittance matrix to use in the test, 4-port
Ys = [  1.701971941826720e-002 - j*5.189648664583783e-002   , -1.047483788864935e-002 + j*2.831281704766298e-002 , ...
      -1.066319031145052e-002 + j*2.849441682290975e-002   , 1.571620548756581e-002  - j*5.662411979813149e-002 ; ...
      -1.047483788864936e-002 + j*2.831281704766301e-002   , 1.701971941826725e-002  - j*5.189648664583781e-002 , ...
 		1.571620548756586e-002 - j*5.662411979813147e-002    , -1.066319031145052e-002  + j*2.849441682290979e-002 ; ...
      -1.066319031145052e-002 +  j*2.849441682290977e-002 , 1.571620548756585e-002 - j*5.662411979813144e-002 , ...
		 1.701971941826726e-002 - j*5.189648664583780e-002   , -1.047483788864935e-002 + j*2.831281704766302e-002 ; ...
       1.571620548756580e-002  - j*5.662411979813149e-002  , -1.066319031145052e-002 + j*2.849441682290975e-002 , ...
   	-1.047483788864935e-002 + j*2.831281704766299e-002  , 1.701971941826719e-002  - j*5.189648664583785e-002 ];

%% Here is how Yc (c is for chain) was obtained
%% 
%% branches = [ 1 2 ; 1 3 ; 1 4 ; 1 5 ; 1 6 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ];

%% nb = size( branches, 1 ); % number of the branches
%% np = size( Ys, 1 ) + 1;        % number of the ports

%% W = zeros( nb, np );
%% K = zeros( nb, np );

%% W( 1:np, 1:np ) = eye( np ); % Voltages applied to ports

%% Y = 1e100*eye( nb,nb ); % ports are the ideal voltage sources
%% Y( np+1:end, np+1:end ) = Ys;

%% [ F, V, I ] = msolve( branches, Y, W, K );

%% % Extracted Y-parameters
%% Yc = -I(1:np,1:np);

Yc = [ ...
[  1.701971941826720e-02 - j*5.189648664583783e-02 , ...
  -2.749455730691656e-02 + j*8.020930369350085e-02 , ...
  -1.883524228011597e-04 + j*1.815997752467602e-04 , ...
   2.637939579901632e-02 - j*8.511853662104127e-02 , ...
  -1.571620548756580e-02 + j*5.662411979813149e-02 ]; ...
[ -2.749455730691655e-02 + j*8.020930369350081e-02 , ... 
   5.498911461383316e-02 - j*1.604186073870016e-01 , ... 
  -1.115161507900241e-03 - j*4.909232927540391e-03 , ... 
  -5.275879159803269e-02 + j*1.702370732420825e-01 , ... 
   2.637939579901632e-02 - j*8.511853662104124e-02 ]; ...
[ -1.883524228011701e-04 + j*1.815997752467706e-04 , ... 
  -1.115161507900218e-03 - j*4.909232927540419e-03 , ... 
   2.607027861402801e-03 + j*9.455266304587304e-03 , ... 
  -1.115161507900243e-03 - j*4.909232927540391e-03 , ... 
  -1.883524228011701e-04 + j*1.815997752467602e-04 ]; ...
[  2.637939579901633e-02 - j*8.511853662104124e-02 , ... 
  -5.275879159803271e-02 + j*1.702370732420825e-01 , ... 
  -1.115161507900230e-03 - j*4.909232927540443e-03 , ... 
   5.498911461383315e-02 - j*1.604186073870016e-01 , ... 
  -2.749455730691654e-02 + j*8.020930369350084e-02 ]; ...
[ -1.571620548756581e-02 + j*5.662411979813149e-02 , ... 
   2.637939579901633e-02 - j*8.511853662104128e-02 , ... 
  -1.883524228011701e-04 + j*1.815997752467706e-04 , ... 
  -2.749455730691654e-02 + j*8.020930369350086e-02 , ... 
   1.701971941826719e-02 - j*5.189648664583785e-02 ] ];

Yt = chainy( Yc );

assertEquals( Ys, Yt, 1.0e-14 );
