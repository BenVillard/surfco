function [varargout] = meshQuality( M , varargin )
% 
% Per cell statistic
% Valid options are:
% for celltype == 3 (polyline)
%     length, l
% for celltype == 5 (triangle mesh)
%     lengths, l
%     minlength, minl
%     maxlength, maxl
%     angles, g
%     minangle, ming
%     maxangle, maxg
%     area, a
%     normal, n
%     heights, h
%     minheight, minh
%     maxheight, maxh
%     inradius, r, ir
%     circumradius, cr
%     aspectratio, ar
%     aspectfrobenius, af
%     edgeratio, er
%     radiusratio, rr
%     relativesize, s
% for celltype == 10 (tetrahedra mesh)
%     lengths
% 
% 
% 

  if size( M.tri,2) > 4
    error('not implemented for cells larger than tetrahedra.');
  end
  
  unitsDEGREE = false;
  
  
  P1 = []; P2 = []; P3 = []; P4 = [];
  %precomputed nodes coordinates
  if size( M.tri,2) > 0
    P1 = M.xyz( M.tri(:,1) ,:); P1(:,end+1:3) = 0;
  end
  if size( M.tri,2) > 1
    P2 = M.xyz( M.tri(:,2) ,:); P2(:,end+1:3) = 0;
  end
  if size( M.tri,2) > 2
    P3 = M.xyz( M.tri(:,3) ,:); P3(:,end+1:3) = 0;
  end
  if size( M.tri,2) > 3
    P4 = M.xyz( M.tri(:,4) ,:); P4(:,end+1:3) = 0;
  end


  fro = @(x) sqrt( sum( x.^2 ,2) );
  nor = @(x) bsxfun( @rdivide , x , sqrt( sum( x.^2 ,2) ) );
  cross = @(a,b)[ a(:,2).*b(:,3) - a(:,3).*b(:,2) ,...
                  a(:,3).*b(:,1) - a(:,1).*b(:,3) ,...
                  a(:,1).*b(:,2) - a(:,2).*b(:,1) ];
  
  L1 = []; L2 = []; L3 = []; L4 = []; L5 = []; L6 = [];
  function get_Ls, if ~isempty( L1 ), return; end
    if size( M.tri,2) > 1
      L1 = P2 - P1;
    end
    if size( M.tri,2) > 2
      L2 = P3 - P2;
      L3 = P1 - P3;
    end
    if size( M.tri,2) > 3
      L4 = P4 - P1;
      L5 = P4 - P2;
      L6 = P4 - P3;
    end
  end

  E1 = []; E2 = []; E3 = []; E4 = []; E5 = []; E6 = [];
  function get_Es, if ~isempty( E1 ), return; end
    get_Ls;

    E1 = nor( L1 );
    E2 = nor( L2 );
    E3 = nor( L3 );
    E4 = nor( L4 );
    E5 = nor( L5 );
    E6 = nor( L6 );
  end

  l1 = []; l2 = []; l3 = []; l4 = []; l5 = []; l6 = []; lm = []; lM = [];
  function get_ls, if ~isempty( l1 ), return; end
    get_Ls;
    
    if size( M.tri,2) > 1
      l1 = fro( L1 );
      lm = l1;
      lM = l1;
    end
    if size( M.tri,2) > 2
      l2 = fro( L2 );
      l3 = fro( L3 );
      lm = min( lm , min( l2 , l3 ) );
      lM = max( lM , max( l2 , l3 ) );
    end
    if size( M.tri,2) > 3
      l4 = fro( L4 );
      l5 = fro( L5 );
      l6 = fro( L6 );
      lm = min( min( lm , l4 ) , min( l5 , l6 ) );
      lM = max( max( lM , l4 ) , max( l5 , l6 ) );
    end
  end

  A1 = []; A2 = []; A3 = []; A4 = [];
  function get_As, if ~isempty( A1 ), return; end
    get_Ls;
    
    if size( M.tri,2) > 2
      A1 = cross( L3 , L1 );
    end
    if size( M.tri,2) > 3
      A1 = -A1;
      A2 = cross( L1 , L5 );
      A3 = cross( L2 , L6 );
      A4 = cross( L3 , L4 );
    end
  end
  
  N1 = []; N2 = []; N3 = []; N4 = [];
  function get_Ns, if ~isempty( N1 ), return; end
    get_As;
    
  	N1 = nor( A1 );
    N2 = nor( A2 );
    N3 = nor( A3 );
    N4 = nor( A4 );
  end
  
  a1 = []; a2 = []; a3 = []; a4 = [];
  function get_as, if ~isempty( a1 ), return; end
    get_As;
    
  	a1 = fro( A1 )/2;
    a2 = fro( A2 )/2;
    a3 = fro( A3 )/2;
    a4 = fro( A4 )/2;
  end

  g1 = []; g2 = []; g3 = [];
  function get_gs, if ~isempty( g1 ), return; end
    get_Es;
  
    if unitsDEGREE
      g1 = 2 * atan2d( fro( E1+E2 ) , fro( E1-E2 ) );
      g2 = 2 * atan2d( fro( E1+E3 ) , fro( E3-E1 ) );
      g3 = 2 * atan2d( fro( E2+E3 ) , fro( E3-E2 ) );
    else
      g1 = 2 * atan2( fro( E1+E2 ) , fro( E1-E2 ) );
      g2 = 2 * atan2( fro( E1+E3 ) , fro( E3-E1 ) );
      g3 = 2 * atan2( fro( E2+E3 ) , fro( E3-E2 ) );
    end
  end
  
  V = [];
  function get_Vs, if ~isempty( V ), return; end
    get_As; get_Ls;
    
    V = dot( A1 , L4 ,2);
  end


  varargout = cell(1,numel(varargin));
  M.celltype = meshCelltype( M );
  for v = 1:numel( varargin )
    if ~ischar( varargin{v} ), error('quality property must be a string.'); end
    switch sprintf( '%s.%d' , lower(varargin{v}) , M.celltype )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'length.3','l.3'}
  get_ls;
  
  varargout{v} = l1;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case {'lengths.5','l.5'}
  get_ls;
  
  varargout{v} = [ l1 , l2 , l3 ];
case {'minlength.5','minl.5'}
  get_ls;
  
  varargout{v} = lm;
case {'maxlength.5','maxl.5'}
  get_ls;
  
  varargout{v} = lM;
case {'angles.5','g.5'}
  get_gs;
  
  varargout{v} = [ g1 , g2 , g3 ];
case {'minangle.5','ming.5'}
  get_gs;
  
  varargout{v} = min( [ g1 , g2 , g3 ] , [] , 2 );
case {'maxangle.5','maxg.5'}
  get_gs;
  
  varargout{v} = max( [ g1 , g2 , g3 ] , [] , 2 );
case {'area.5','a.5'}
  get_as;
  
  varargout{v} = a1;
case {'normal.5','n.5'}
  get_Ns;
  
  varargout{v} = N1;
case {'heights.5','h.5'}
  get_as; get_ls;
  
  varargout{v} = bsxfun( @rdivide , a1 , [l1 l2 l3] )*2;
case {'minheight.5','minh.5'}
  get_as; get_ls;
  
  varargout{v} = min( bsxfun( @rdivide , a1 , [l1 l2 l3] ) , [] , 2 )*2;
case {'maxheight.5','maxh.5'}
  get_as; get_ls;
  
  varargout{v} = max( bsxfun( @rdivide , a1 , [l1 l2 l3] ) , [] , 2 )*2;
case {'inradius.5','r.5','ir.5','in.5'}
  get_as; get_ls;

  varargout{v} = 2 * a1 ./ ( l1 + l2 + l3 );
case {'circumradius.5','cr.5','circumr.5'}
  get_as; get_ls;

  varargout{v} = ( l1 .* l2 .* l3 ) ./ ( 4 * a1 );
case {'aspectratio.5','ar.5'}
  get_as; get_ls;

  varargout{v} = lM .* ( l1 + l2 + l3 ) ./ ( 4 * sqrt(3) * a1 );
case {'aspectfrobenius.5','af.5'}
  get_as; get_ls;

  varargout{v} = ( l1.^2 + l2.^2 + l3.^2 ) ./ ( 4 * sqrt(3) * a1 );
case {'edgeratio.5','er.5'}
  get_ls;

  varargout{v} = lM ./ lm;
case {'radiusratio.5','rr.5'}
  get_ls; get_as;

  varargout{v} = ( l1 .* l2 .* l3 ) .* ( l1 + l2 + l3 ) ./ ( 16 * a1.^2 );
case {'relativesize.5','s.5','rs.5'}
  get_as;
  
  r = a1 ./ mean( a1 );
  varargout{v} = min( r , 1./r );
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case {'lengths.10'}
  get_ls;
  
  varargout{v} = [ l1 , l2 , l3 , l4 , l5 , l6 ];
case {'areas.10'}
  get_as;
  
  varargout{v} = [ a1 , a2 , a3 , a4 ];
case {'volume.10'}
  get_Vs;
  
  varargout{v} = abs( V/6 );

case {'signedvolume.10'}
  get_As; get_Ls;
  
  varargout{v} = - dot( A1 , L4 ,2)/6;
  
case {'orientation.10'}
  get_As; get_Ls;

  varargout{v} = -sign( dot( A1 , L4 ,2) );

  
  

        case {'inradius.10'}
          V = dot( cross( L3 , L1 ) , L4 ,2)/6;
          A = ( fro( cross( L3 , L1 ) ) + fro( cross( L4 , L1 ) ) + fro( cross( L5 , L2 ) ) + fro( cross( L4 , L3 ) ) )/2;

          varargout{v} = 3 * V ./ A;
          
        case {'circumradius.10'}
          V = dot( cross( L3 , L1 ) , L4 ,2)/6;
          varargout{v} = fro( bsxfun( @times , l4.^2 , cross( L3 , L1 ) ) + ...
                              bsxfun( @times , l3.^2 , cross( L4 , L1 ) ) + ...
                              bsxfun( @times , l1.^2 , cross( L4 , L3 ) ) ) ./ ( 12 * V );


        case {'aspectbeta.10'}

        case {'aspectgamma.10'}%OK
          V = dot( cross( L3 , L1 ) , L4 ,2)/6;

          varargout{v} = sqrt( ( l1.^2 + l2.^2 + l3.^2 + l4.^2 + l5.^2 + l6.^2 )/6 ).^3 * sqrt(2) ./ ( 12 * V );

        case {'aspectfrobenius.10'}

        case {'aspectratio.10'}%OK
          V = dot( cross( L3 , L1 ) , L4 ,2)/6;
          A = ( fro( cross( L3 , L1 ) ) + fro( cross( L4 , L1 ) ) + fro( cross( L5 , L2 ) ) + fro( cross( L4 , L3 ) ) )/2;
          r = 3 * V ./ A;

          varargout{v} = lM ./ ( 2*sqrt(6)*r );

        case {'collapseratio.10'}
          H = @(A,b) abs( dot(A,b,2) ./ fro( A ) );
          h0 = H( cross( L2 , L5 ) , L1 );
          h1 = H( cross( L3 , L4 ) , L1 );
          h2 = H( cross( L1 , L4 ) , L2 );
          h3 = H( cross( L1 , L2 ) , L4 );
          
          varargout{v} = min( [ h0 ./ max( [ l2 , l5 , l6 ] , [] , 2 ) ,...
                                h1 ./ max( [ l3 , l4 , l6 ] , [] , 2 ) ,...
                                h2 ./ max( [ l1 , l4 , l5 ] , [] , 2 ) ,...
                                h3 ./ max( [ l1 , l2 , l3 ] , [] , 2 ) ] , [] , 2 );
          
        case {'condition.10'}

        case {'distortion.10'}

        case {'jacobian.10'} %OK
          varargout{v} = dot( cross( L3 , L1 ) , L4 ,2);

        case {'dihedralangles.10'}
%           ang = @(b1,b2,b3)atan2( ...
%             dot( cross( cross(b1,b2) , cross(b2,b3) ) , nor( b2 ) , 2 ) ,...
%             dot( cross( b1 , b2 ) , cross( b2 , b3 ) ,2) );
%           
%           a0 = ang( L0 , L1 , L4 );
%           a1 = ang( L1 , L0 , L4 );
%           a2 = ang( L2 , L0 , L3 );
%           a3 = ang( L3 , L0 , L2 );
%           a4 = ang( L4 , L0 , L1 );
%           a5 = ang( L5 , L1 , L2 );
          ang = @(na,nb) acos( - dot( na , nb ,2 ) ./ ( fro(na) .* fro(nb) ) );
          n0 = cross( L2 , L5 );
          n1 = cross( L3 , L4 );
          n2 = cross( L1 , L4 );
          n3 = cross( L2 , L1 );
          
          a1 = ang( n3  , n2  );
          a2 = ang( n3  , n0  );
          a3 = ang( n3  , n1  );
          a4 = ang( n2  , n1  );
          a4 = ang( n0  , n2  );
          a5 = ang( n0  , n1  );
          
          varargout{v} = max( [ a1 , a2 , a3 , a4 , a4 , a5 ] , [] , 2 ) * 180/pi;

% % % % % case {'lmin.3','minedge.3','lmin.5','minedge.5','lmin.10','minedge.10'}
% % % % %   varargout{v} = lmin;
% % % % % case {'lmax.3','maxedge.3','lmax.5','maxedge.5','lmax.10','maxedge.10'}
% % % % %   varargout{v} = lmax;
% % % % % case {'edgeratio.5','edgeratio.10'}
% % % % %   varargout{v} = lmax ./ lmin;
% % % % %         case {'radiusratio.10'} %OK
% % % % %           V = dot( cross( L2 , L0 ) , L3 ,2)/6;
% % % % %           A = ( fro( cross( L2 , L0 ) ) + fro( cross( L3 , L0 ) ) + fro( cross( L4 , L1 ) ) + fro( cross( L3 , L2 ) ) )/2;
% % % % %           r = 3 * V ./ A;
% % % % %           R = fro(  bsxfun( @times , l3.^2 , cross( L2 , L0 ) ) + ...
% % % % %                     bsxfun( @times , l2.^2 , cross( L3 , L0 ) ) + ...
% % % % %                     bsxfun( @times , l0.^2 , cross( L3 , L2 ) ) ) ./ ( 12 * V );
% % % % % 
% % % % %           varargout{v} = R ./ ( 3*r );

        case {'relativesizesquared.10'}

        case {'scaledjacobian.10'}

        case {'shape.10'}

        case {'shapeandsize.10'}
          
        case {'aspectdelta.10'}


          
          
      otherwise, error('invalid property "%s" for celltype %d (or not implemented yet! :P).',varargin{v},M.celltype);
    end
    
  end
          





end
