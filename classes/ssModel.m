%% ssModel - Class for space state models
%
% This class is an extension of the "ss" class. It is designed as an extended
% vay to define state space models,suitable for working with systems
% obtained from linealization of a non-lienar system arround an operating point.
%
% Let us consider a discrete time-invariant state space model
%
%       x+ = A*x + B*in                                     (1a)
%       y  = C*x + D*in                                     (1b)
%
% where "in" are all the input signals of the system and "y" all the output signals.
% The input signals can be divided in three groups: the control inputs "u", the
% unmeasured disturbances "du" and the measured disturbances "dm".
% The output signals can also be divided in three groups: the constrained outputs "yc",
% the tracked outputs "yr" and the measured outputs "ym".
% We assume that this model is describing the dynamics of a non-linear system arround
% an operating point given by the state X0, the inputs IN0 and the outputs Y0.
% 
% This class is used to represent systems described by (1). It is given the matrices
% A, B, C and D, as well as vectors, u, du, dm, yc, yr and ym indicating which columns/rows
% of the matrices B, C and D correspond to each one of the signals. Additionally, it also
% requires the operating point (X0, IN0, Y0).
% It then describes model (1) as the following discrete time-invariant state space model
%
%       x+ = A*x + Bu*u + Bdu*du + Bdm*dm                   (2a)
%       yc = Cc*x + Dcu*u + Dcdu*du + Dcdm*dm               (2b)
%       yr = Cr*x + Dru*u + Drdu*du + Drdm*dm               (2c)
%       ym = Cm*x + Dmu*u + Dmdu*du + Dmdm*dm,              (2d)
%
% We consider that the system can be subject to the following constraints,
% 
%       LBx <=  x  <= UBx
%       LBu <=  u  <= UBu
%       LBy <= y_c <= UBy
%
% The class distinguished between "engineering units", which are expressed in the units
% of the non-linear model, and "incremental units", which are expressed in incremental
% units with respect to the operating point.
%
% Finally, the class allows a normalization of the linear model, where the state vector "x",
% the input vector "in" and the output vector "y" are subject to the following normalization
% via the scaling vectors xN, inN and yN:
% 
%       X = (1./xN).*x + X0, IN = (1./inN).*in + IN0, Y = (1./yN).*y + Y0, (3)
%
% where X, IN and Y are the state, inputs and outputs in engineering units.
% 
% Constructor: sys = ssModel(A, B, C, D, Ts, uu, du, dm, yr, yc, ym, LBxEng, UBxEng, LBuEng, UBuEng, LByEng, UByEng, X0, IN0, Y0, xN, inN, yN)
%    - (A, B, C, D) are required
%   - Ts is the sampling time of the system. If it is not provided, then ssModel
%     assumes that A, B, C and D describes a discrete model of unknown sample time.
%     For a continuous-time model: sys = ssModel(A, B, C, D, 0)
%     In any case, transformation between continuous-time and discrete-time
%     models can be done using the functions c2d and d2c, as with any ss instance.
%   - The constructor can also be given the following optional arguments:
%     (uu, du, dm, yr, yc, ym, LBxEng, UBxEng, LBuEng, UBuEng, LByEng, UByEng, X0, IN0, Y0, xN, inN, yN)
%     These properties of this class can be updated after it is constructed, but doing so can have unnexpected
%     consequences. In order to avoid them, it is important that they be updated in the following order:
%       - Operating point: X0, IN0, Y0
%       - Bounds in engineering units: LBxEng, UBxEng, LBuEng, UBuEng, LByEng, UByEng
%       - Scaling vectors: Nx, Nu, Ndu, Ndm, Nyc, Nyr, Nym
%         We note that the scaling vectors xN, inN and yN can only be updated by the constructor.
%     Avoid changing the variables for the bounds in incremental units (LBx, UBx, LBu, UBu, LBy and UBy) if possible,
%     unless you realli know what you are doing, since it can have unexpected consecuences.
%
% Important properties:
%   - X0, IN0, Y0, u0, du0, dm0, yc0, yr0, ym0: Operating point
%   - LBx, UBx, LBu, UBu, LBy, UBy: Constraints in incremental units
%   - LBxEng, UBxEng, LBuEng, UBuEng, LByEng, UByEng: Constraints in engineering units
%   - uu: Columns of B that belong to the system inputs
%   - du: Columns of B that belong to the system unmeasured disturbances
%   - dm: Columns of B that belong to the system measured disturbances
%   - yr: Rows of C that belong to the tracked outputs
%   - yc: Rows of C that belong to the constrained outputs
%   - ym: Rows of C that belong to the measured outputs
%   - A, B, C, D: Matrices of model (1)
%   - A, Bu, Bdu, Bdm, Cc, Dcu, Dcm, etc: Matrices of model (2)
%   - n_x, n_u, n_du, n_dm, n_yc, n_yr, n_ym: Dimension of each signal
%   - Nx, Nu, Ndu, Ndm, Nyc, Nyr, Nym: Scaling vectors for each signal
%
% All the above properties can be set, accesed and changed after construction.
%
% Methods:
%   - x = x2Inc(X): Returnes the state x in incremental units corresponding
%     to the state X, given in engineering units. This takes into account the 
%     operating point and the normalization vectors as in (3).
%   - X = x2Eng(x): Same as x2Inc, but in the other direction.
%   - Similar methods exist for all the other appropiate properties of the class.
%
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

classdef ssModel < ss
    %% PROPERTIES
    properties
        LBxEng % Lower bounds for the system state. In engineering units
        UBxEng % Upper bounds for the system state. In engineering units
        LBuEng % Lower bounds for the system input. In engineering units
        UBuEng % Upper bounds for the system input. In engineering units
        LByEng % Lower bounds for the system constrained output. In engineering units
        UByEng % Upper bounds for the system constrained output. In engineering units
        X0 % Operating point of system state
        IN0 % Operating point of system inputs
        Y0 % Operating point of system outputs
    end
    properties (Hidden=true)
        % These scaling vectors are hidden because they should not be changed except by using method "normalize"
        xN % Scaling vector for the system states
        inN % Scaling vector for the system inputs and disturbances
        yN % Scaling vector for the system outputs
    end
    properties (Dependent, Hidden=false)
        uu  % Columns of B that belong to the system inputs
        du %  Columns of B that belong to the system unmeasured disturbances
        dm %  Columns of B that belong to the system measured disturbances
        yr % Rows of C that belong to the tracked outputs
        yc % Rows of C that belong to the constrained outputs
        ym % Rows of C that belong to the measured outputs
        u0 % Operating point for the control action
        dm0 % Operating point for the measured disturbance
        du0 % Operating point for the unmeasured disturbance
        yr0 % Operating point for the tracked output
        yc0 % Operating point for the constrained output
        ym0 % Operating point for the measured output
        Bu % B matrix for the system inputs
        Bdu % B matrix for the system unmeasured disturbances
        Bdm % B matrix for the system measured disturbances
        Cr % C matric of the model for the tracked outputs
        Cc % C matrix of the model for the constrained outputs
        Cm % C matrix of the model for the measured outputs
        Dru % D matric of the model for the tracked outputs from system inputs
        Drdu % D matric of the model for the tracked outputs from unmeasured disturbances
        Drdm % D matric of the model for the tracked outputs from measured disturbances
        Dcu % D matrix of the model for the constrained outputs from system inputs
        Dcdu % D matric of the model for the constrained outputs from unmeasured disturbances
        Dcdm % D matric of the model for the constrained outputs from measured disturbances
        Dmu % D matrix of the model for the measured outputs from system inputs
        Dmdu % D matric of the model for the measured outputs from unmeasured disturbances
        Dmdm % D matric of the model for the measured outputs from measured disturbances
        n_x % Number of system states
        n_in % Number of system inputs
        n_y % Number of system outputs
        n_u % Number of system control actions
        n_du % Number of unmeasured disturbances
        n_dm % Number of measured disturbances
        n_yr % Number of tracked outputs
        n_yc % Number of constrained outputs
        n_ym % Number of measured outputs
        Nx % Scaling vector for the system states
        Nin % Scaling vector for system inputs (all of them)
        Ny % Scaling vector for the system outputs (all of them)
        Nu % Scaling vector for the inputs (control actions)
        Ndu % Scaling vector for unmeasured disturbances
        Ndm % Scaling vector for measured disturbances
        Nyr % Scaling vector for tracked outputs
        Nyc % Scaling vector for constrained outputs
        Nym % Scaling vector for measured outputs
        LBx % Lower bounds for the system state. In incremental and scaled value
        UBx % Upper bounds for the system state. In incremental and scaled value
        LBu % Lower bounds for the system input. In incremental and scaled value
        UBu % Upper bounds for the system input. In incremental and scaled value
        LBy % Lower bounds for the system constrained output. In incremental value
        UBy % Upper bounds for the system constrained output. In incremental value
    end
    properties (Dependent, Hidden=true)
        isCtrb % Boolean that indicates if the system is controllable
        isObsv % Boolean that indicates if the system is observable
    end
    
    methods
    %% CONSTRUCTOR
    function self = ssModel(A, B, C, D, varargin) % (A, B, C, D, Ts, uu, du, dm, yr, yc, ym, LBxEng, UBxEng, LBuEng, UBuEng, LByEng, UByEng, X0, IN0, Y0, xN, inN, yN)
        % Defaults
        nx = size(A, 1);
        nu = size(B, 2);
        ny = size(C, 1);
        def_Ts = []; % Note: This defaults the constructor to creating distrete-time models with unknown sample time. The ss constructor will change this to Ts = -1.
        def_uu = 1:nu;
        def_du = [];
        def_dm = [];
        def_yr = 1:ny;
        def_yc = 1:ny;
        def_ym = 1:ny;
        def_LBxEng = -inf*ones(nx, 1);
        def_UBxEng = inf*ones(nx, 1);
        def_LBuEng = -inf*ones(nu, 1);
        def_UBuEng = inf*ones(nu, 1);
        def_LByEng = [];
        def_UByEng = [];
        def_X0 = zeros(nx, 1);
        def_IN0 = zeros(nu, 1);
        def_Y0 = zeros(ny, 1);
        def_xN = ones(nx, 1);
        def_uN = ones(nu, 1);
        def_yN = ones(ny, 1);
        % Parser
        p = inputParser;
        p.CaseSensitive = false;
        p.FunctionName = 'ssModel_Constructor';
        % Required
        addRequired(p, 'A', @(x) isnumeric(x));
        addRequired(p, 'B', @(x) isnumeric(x));
        addRequired(p, 'C', @(x) isnumeric(x));
        addRequired(p, 'D', @(x) isnumeric(x));
        % Optional
        addOptional(p, 'Ts', def_Ts, @(x) (isnumeric(x) && isscalar(x))||isempty(x));
        addOptional(p, 'uu', def_uu, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'du', def_du, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'dm', def_dm, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'yr', def_yr, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'yc', def_yc, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'ym', def_ym, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'LBxEng', def_LBxEng, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nx) || isempty(x));
        addOptional(p, 'UBxEng', def_UBxEng, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nx) || isempty(x));
        addOptional(p, 'LBuEng', def_LBuEng, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nu) || isempty(x));
        addOptional(p, 'UBuEng', def_UBuEng, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nu) || isempty(x));
        addOptional(p, 'LByEng', def_LByEng, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'UByEng', def_UByEng, @(x) (isnumeric(x) && (min(size(x))==1)) || isempty(x));
        addOptional(p, 'X0', def_X0, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nx) || isempty(x));
        addOptional(p, 'IN0', def_IN0, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nu) || isempty(x));
        addOptional(p, 'Y0', def_Y0, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==ny) || isempty(x));
        addOptional(p, 'xN', def_xN, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nx) || isempty(x));
        addOptional(p, 'inN', def_uN, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==nu) || isempty(x));
        addOptional(p, 'yN', def_yN, @(x) (isnumeric(x) && (min(size(x))==1) && length(x)==ny) || isempty(x));
        % Parse
        parse(p, A, B, C, D, varargin{:})
        % Rename
        uu = p.Results.uu;
        du = p.Results.du;
        dm = p.Results.dm;
        yr = p.Results.yr;
        yc = p.Results.yc;
        ym = p.Results.ym;
        LBxEng = p.Results.LBxEng;
        UBxEng = p.Results.UBxEng;
        LBuEng = p.Results.LBuEng;
        UBuEng = p.Results.UBuEng;
        LByEng = p.Results.LByEng;
        UByEng = p.Results.UByEng;
        X0 = p.Results.X0;
        IN0 = p.Results.IN0;
        Y0 = p.Results.Y0;
        xN = p.Results.xN;
        inN = p.Results.inN;
        yN = p.Results.yN;
        % Defauls
        if isempty(yr); yr = def_yr; end
        if isempty(yc); yc = def_yc; end
        if isempty(ym); ym = def_ym; end
        if isempty(LBxEng); LBxEng = def_LBxEng; end
        if isempty(UBxEng); UBxEng = def_UBxEng; end
        if isempty(LBuEng); LBuEng = def_LBuEng; end
        if isempty(UBuEng); UBuEng = def_UBuEng; end
        if isempty(X0); X0 = def_X0; end
        if isempty(IN0); IN0 = def_IN0; end
        if isempty(Y0); Y0 = def_Y0; end
        % Make vectors column or row vectors
        if size(uu,1)>1; uu=uu'; end % Make u a row vector
        if size(du,1)>1; du=du'; end % Make du a row vector
        if size(dm,1)>1; dm=dm'; end % Make dm a row vector
        if size(yr,1)>1; yr=yr'; end % Make yr a row vector
        if size(yc,1)>1; yc=yc'; end % Make yc a row vector
        if size(ym,1)>1; ym=ym'; end % Make ym a row vector
        if size(X0,2)>1; X0=X0'; end % Make LBx a column vector
        if size(IN0,2)>1; IN0=IN0'; end % Make LBx a column vector
        if size(Y0,2)>1; Y0=Y0'; end % Make LBx a column vector
        if size(LBxEng,2)>1; LBxEng=LBxEng'; end % Make LBx a column vector
        if size(UBxEng,2)>1; UBxEng=UBxEng'; end % Make UBx a column vector
        if size(LBuEng,2)>1; LBuEng=LBuE'; end % Make LBu a column vector
        if size(UBuEng,2)>1; UBuEng=UBuEng'; end % Make UBu a column vector
        if size(LByEng,2)>1; LByEng=LByEng'; end % Make LBy a column vector
        if size(UByEng,2)>1; UByEng=UByEng'; end % Make UBy a column vector
        % Perform checks
        if isempty(p.Results.LByEng)
            LByEng = -inf*ones(length(yc), 1);
        else
            if length(LByEng) ~= length(yc); warning('ssModel:DimensionError:LBy', 'Length of LBy must be equal to length of yc'); end
        end
        if isempty(p.Results.UByEng)
            UByEng = inf*ones(length(yc), 1);
        else
            if length(UByEng) ~= length(yc); warning('ssModel:DimensionError:UBy', 'Length of UBy must be equal to length of yc'); end
        end

        % Construct ss object
        self@ss(p.Results.A, p.Results.B, p.Results.C, p.Results.D, p.Results.Ts)
        self.InputGroup.uu = uu;
        self.InputGroup.du = du;
        self.InputGroup.dm = dm;
        self.OutputGroup.yr = yr;
        self.OutputGroup.yc = yc;
        self.OutputGroup.ym = ym;

        % Set class properties
        self.LBxEng = LBxEng;
        self.UBxEng = UBxEng;
        self.LBuEng = LBuEng;
        self.UBuEng = UBuEng;
        self.LByEng = LByEng;
        self.UByEng = UByEng;
        self.X0  = X0;
        self.IN0 = IN0;
        self.Y0  = Y0;
        self.xN = def_xN;
        self.inN = def_uN;
        self.yN = def_yN;
        % Normalize
        self = normalize(self, xN, inN, yN);

    end
    
    %% SETTERS and GETTERS
    
    % Input groups
    % uu
    function r = get.uu(self); try r = self.InputGroup.uu; catch; r = []; end; end % Getter
    function self = set.uu(self, value) % Setter
        if size(value,1)>1; value=value'; end % Make a row vector
        self.InputGroup.uu = value;
        self.LBu = -Inf*ones(length(value), 1);
        self.UBu = Inf*ones(length(value), 1);
    end
    % du
    function r = get.du(self); try r = self.InputGroup.du; catch; r = []; end; end % Getter
    function self = set.du(self, value) % Setter
        if size(value,1)>1; value=value'; end % Make a row vector
        self.InputGroup.du = value;
    end
    % dm
    function r = get.dm(self); try r = self.InputGroup.dm; catch; r = []; end; end % Getter
    function self = set.dm(self, value) % Setter
        if size(value,1)>1; value=value'; end % Make a row vector
        self.InputGroup.dm = value;
    end
    
    % B matrices
    function r = get.Bu(self); r = self.B(:, self.uu); end % Bu getter
    function r = get.Bdu(self); r = self.B(:, self.du); end % Bdu getter
    function r = get.Bdm(self); r = self.B(:, self.dm); end % Bdm getter
    
    % Output groups
    % yr
    function r = get.yr(self); try r = self.OutputGroup.yr; catch; r = []; end; end % Getter
    function self = set.yr(self, value) % Setter
        if size(value,1)>1; value=value'; end % Make a row vector
        self.OutputGroup.yr = value;
    end
    % yc
    function r = get.yc(self); try r = self.OutputGroup.yc; catch; r = []; end; end % Getter
    function self = set.yc(self, value) % Setter
        if size(value,1)>1; value=value'; end % Make a row vector
        self.OutputGroup.yc = value;
        self.LBy = -Inf*ones(length(value), 1);
        self.UBy = Inf*ones(length(value), 1);
    end
    % ym
    function r = get.ym(self); try r = self.OutputGroup.ym; catch; r = []; end; end % Getter
    function self = set.ym(self, value) % Setter
        if size(value,1)>1; value=value'; end % Make a row vector
        self.OutputGroup.ym = value;
    end
    
    % C matrices
    function r = get.Cr(self); r = self.C(self.yr, :); end % Getter
    function r = get.Cc(self); r = self.C(self.yc, :); end % Getter
    function r = get.Cm(self); r = self.C(self.ym, :); end % Getter

    % D matrices
    % Dr
    function r = get.Dru(self); r = self.D(self.yr, self.uu); end % Dru getter
    function r = get.Drdu(self); r = self.D(self.yr, self.du); end % Drdu getter
    function r = get.Drdm(self); r = self.D(self.yr, self.dm); end % Drdm getter
    % Dc
    function r = get.Dcu(self); r = self.D(self.yc, self.uu); end % Dcu getter
    function r = get.Dcdu(self); r = self.D(self.yc, self.du); end % Dcdu getter
    function r = get.Dcdm(self); r = self.D(self.yc, self.dm); end % Dcdm getter
    % Dm
    function r = get.Dmu(self); r = self.D(self.ym, self.uu); end % Dmu getter
    function r = get.Dmdu(self); r = self.D(self.ym, self.du); end % Dmdu getter
    function r = get.Dmdm(self); r = self.D(self.ym, self.dm); end % Dmdm getter
    
    % Operating point
    % u0
    function r = get.u0(self); r = self.IN0(self.uu); end % Setter
    function self = set.u0(self, value) % Getter
        if size(value,2)>1; value=value'; end % Make column vector
        self.IN0(self.uu) = value;
    end
    % dm0
    function r = get.dm0(self); r = self.IN0(self.dm); end % Setter
    function self = set.dm0(self, value) % Getter
        if size(value,2)>1; value=value'; end % Make column vector
        self.IN0(self.dm) = value;
    end
    % du0
    function r = get.du0(self); r = self.IN0(self.du); end % Setter
    function self = set.du0(self, value) % Getter
        if size(value,2)>1; value=value'; end % Make column vector
        self.IN0(self.du) = value;
    end
    % yr0
    function r = get.yr0(self); r = self.Y0(self.yr); end % Setter
    function self = set.yr0(self, value) % Getter
        if size(value,2)>1; value=value'; end % Make column vector
        self.Y0(self.yr) = value;
    end
    % yc0
    function r = get.yc0(self); r = self.Y0(self.yc); end % Setter
    function self = set.yc0(self, value) % Getter
        if size(value,2)>1; value=value'; end % Make column vector
        self.Y0(self.yc) = value;
    end
    % ym0
    function r = get.ym0(self); r = self.Y0(self.ym); end % Setter
    function self = set.ym0(self, value) % Getter
        if size(value,2)>1; value=value'; end % Make column vector
        self.Y0(self.ym) = value;
    end
    
    % Mayor scaling vectors
    % Nx
    function r = get.Nx(self); r = self.xN; end % Getter
    function self = set.Nx(self, value); self = normalize(self, value, self.inN, self.yN); end % Setter
    % Nu
    function r = get.Nin(self); r = self.inN; end % Getter
    function self = set.Nin(self, value); self = normalize(self, self.xN, value, self.yN); end % Setter
    % Ny
    function r = get.Ny(self); r = self.yN; end % Getter
    function self = set.Ny(self, value); self = normalize(self, self.xN, self.inN, value); end % Setter
    
    % Minor scaling vectors
    % Nuu
    function r = get.Nu(self); r = self.inN(self.uu); end % Getter
    function self = set.Nu(self, value) % Setter
        aux_uN = self.inN;
        aux_uN(self.uu) = value;
        self = normalize(self, self.xN, aux_uN, self.yN);
    end
    % Ndu
    function r = get.Ndu(self); r = self.inN(self.du); end % Getter
    function self = set.Ndu(self, value) % Setter
        aux_uN = self.inN;
        aux_uN(self.du) = value;
        self = normalize(self, self.xN, aux_uN, self.yN);
    end
    % Ndm
    function r = get.Ndm(self); r = self.inN(self.dm); end % Getter
    function self = set.Ndm(self, value) % Setter
        aux_uN = self.inN;
        aux_uN(self.dm) = value;
        self = normalize(self, self.xN, aux_uN, self.yN);
    end
    % Nyr
    function r = get.Nyr(self); r = self.yN(self.yr); end % Getter
    function self = set.Nyr(self, value) % Setter
        aux_yN = self.yN;
        aux_yN(self.yr) = value;
        self = normalize(self, self.xN, self.inN, aux_yN);
    end
    % Nyc
    function r = get.Nyc(self); r = self.yN(self.yc); end % Getter
    function self = set.Nyc(self, value) % Setter
        aux_yN = self.yN;
        aux_yN(self.yc) = value;
        self = normalize(self, self.xN, self.inN, aux_yN);
    end
    % Nym
    function r = get.Nym(self); r = self.yN(self.ym); end % Getter
    function self = set.Nym(self, value) % Setter
        aux_yN = self.yN;
        aux_yN(self.ym) = value;
        self = normalize(self, self.xN, self.inN, aux_yN);
    end
    
    % Dimensions
    function n_x = get.n_x(self); n_x = size(self.A, 1); end
    function n_in = get.n_in(self); n_in = size(self.B, 2); end
    function n_y = get.n_y(self); n_y = size(self.C, 1); end
    function n_u = get.n_u(self); try n_u = length(self.InputGroup.uu); catch; n_u = 0; end; end
    function n_du = get.n_du(self); try n_du = length(self.InputGroup.du); catch; n_du = 0; end; end
    function n_dm = get.n_dm(self); try n_dm = length(self.InputGroup.dm); catch; n_dm = 0; end; end
    function n_yr = get.n_yr(self); try n_yr = length(self.OutputGroup.yr); catch; n_yr = 0; end; end
    function n_yc = get.n_yc(self); try n_yc = length(self.OutputGroup.yc); catch; n_yc = 0; end; end
    function n_ym = get.n_ym(self); try n_ym = length(self.OutputGroup.ym); catch; n_ym = 0; end; end
    
    % Controllable and observable
    function r = get.isCtrb(self)
        % Returns true is the system is controllable
        r = rank(ctrb(self)) == self.n_x;
    end
    function r = get.isObsv(self)
        % Returs true if the system is observable
        r = rank(Obsv(self)) == self.n_x;
    end
    
    % Bounds
    % LBx 
    function r = get.LBx(self); r = x2Inc(self, self.LBxEng); end % Getter
    function self = set.LBx(self, value) % Setter
        if size(value,2)>1; value=value'; end % Make column vector
        self.LBxEng = x2Eng(self, value);
    end
    % UBx 
    function r = get.UBx(self); r = x2Inc(self, self.UBxEng); end % Getter
    function self = set.UBx(self, value) % Setter
        if size(value,2)>1; value=value'; end % Make column vector
        self.UBxEng = x2Eng(self, value);
    end
    % LBu
    function r = get.LBu(self); r = u2Inc(self, self.LBuEng); end % Getter
    function self = set.LBu(self, value) % Setter
        if size(value,2)>1; value=value'; end % Make column vector
        self.LBuEng = u2Eng(self, value);
    end
    % UBu
    function r = get.UBu(self); r = u2Inc(self, self.UBuEng); end % Getter
    function self = set.UBu(self, value) % Setter
        if size(value,2)>1; value=value'; end % Make column vector
        self.UBuEng = u2Eng(self, value);
    end
    % LBy
    function r = get.LBy(self); r = yc2Inc(self, self.LByEng); end % Getter
    function self = set.LBy(self, value) % Setter
        if size(value,2)>1; value=value'; end % Make column vector
        self.LByEng = yc2Eng(self, value);
    end
    % UBy in engineering units
    function r = get.UBy(self); r = yc2Inc(self, self.UByEng); end % Getter
    function self = set.UBy(self, value) % Setter
        if size(value,2)>1; value=value'; end % Make column vector
        self.UByEng = yc2Eng(self, value);
    end
        
    %% METHODS    
    function O = Obsv(self)
        % Returns the observability matrix for the measured outputs (i.e. uses Cm instead of for all C).
        O = obsv(self.A, self.Cm); 
    end
    
    function self = normalize(self, Nx, Nu, Ny)
        % Normalizes the system model for new values of vectors Nx, Nu and Ny. Updates the properties accordingly
        % Calculate and update new system model matrices
        self.A = diag(Nx)*diag(1./self.xN)*self.A*diag(self.xN)*diag(1./Nx);
        self.B = diag(Nx)*diag(1./self.xN)*self.B*diag(self.inN)*diag(1./Nu);
        self.C = diag(Ny)*diag(1./self.yN)*self.C*diag(self.xN)*diag(1./Nx);
        self.D = diag(Ny)*diag(1./self.yN)*self.D*diag(self.inN)*diag(1./Nu);
        % Update values of the scaling vectors
        self.xN = Nx;
        self.inN = Nu;
        self.yN = Ny;
    end
    
    % Transforms a scaled and incremental x into one in engineering units
    function xEng = x2Eng(self, xInc)
        xEng = zeros(self.n_x, size(xInc, 2));
        for i=1:size(xInc, 2)
            xEng(:,i) = (1./self.Nx).*xInc(:,i) + self.X0;
        end
    end
    % Transforms an x in engineering units into one in scaled and incremental units
    function xInc = x2Inc(self, xEng)
        xInc = zeros(self.n_x, size(xEng, 2));
        for i=1:size(xEng, 2)
            xInc(:,i) = self.Nx.*(xEng(:,i) - self.X0);
        end
    end
    % Transforms a scaled and incremental u (control action) into one in engineering units
    function uEng = u2Eng(self, uInc)
        uEng = zeros(self.n_u, size(uInc, 2));
        for i=1:size(uInc, 2)
            uEng(:,i) = (1./self.Nu).*uInc(:,i) + self.u0;
        end
    end
    % Transforms a u in engineering units into one in scaled and incremental units
    function uInc = u2Inc(self, uEng)
        uInc = zeros(self.n_u, size(uEng, 2));
        for i=1:size(uEng, 2)
            uInc(:,i) = self.Nu.*(uEng(:,i) - self.u0);
        end
    end
    % Transforms a scaled and incremental dm  into one in engineering units
    function dmEng = dm2Eng(self, dmInc)
        dmEng = zeros(self.n_dm, size(dmInc, 2));
        for i=1:size(dmInc, 2)
            dmEng(:,i) = (1./self.Ndm).*dmInc(:,i) + self.dm0;
        end
    end
    % Transforms a dm in engineering units into one in scaled and incremental units
    function dmInc = dm2Inc(self, dmEng)
        dmInc = zeros(self.n_dm, size(dmEng, 2));
        for i=1:size(dmEng, 2)
            dmInc(:,i) = self.Ndm.*(dmEng(:,i) - self.dm0);
        end
    end
    % Transforms a scaled and incremental du  into one in engineering units
    function duEng = du2Eng(self, duInc)
        duEng = zeros(self.n_du, size(duInc, 2));
        for i=1:size(duInc, 2)
            duEng(:,i) = (1./self.Ndu).*duInc(:,i) + self.du0;
        end
    end
    % Transforms a du in engineering units into one in scaled and incremental units
    function duInc = du2Inc(self, duEng)
        duInc = zeros(self.n_du, size(duEng, 2));
        for i=1: size(duEng, 2)
            duInc(:,i) = self.Ndu.*(duEng(:,i) - self.du0);
        end
    end
    % Transforms a scaled and incremental yr into one in engineering units
    function yrEng = yr2Eng(self, yrInc)
        yrEng = zeros(self.n_yr, size(yrInc, 2));
        for i=1:size(yrInc, 2)
            yrEng(:,i) = (1./self.Nyr).*yrInc(:,i) + self.yr0;
        end
    end
    % Transforms an yr in engineering units into one in scaled and incremental units
    function yrInc = yr2Inc(self, yrEng)
        yrInc = zeros(self.n_yr, size(yrEng, 2));
        for i=1:size(yrEng, 2)
            yrInc(:,i) = self.Nyr.*(yrEng(:,i) - self.yr0);
        end
    end
    % Transforms a scaled and incremental yc into one in engineering units
    function ycEng = yc2Eng(self, ycInc)
        ycEng = zeros(self.n_yc, size(ycInc, 2));
        for i=1:size(ycInc, 2)
            ycEng(:,i) = (1./self.Nyc).*ycInc(:,i) + self.yc0;
        end
    end
    % Transforms an yc in engineering units into one in scaled and incremental units
    function ycInc = yc2Inc(self, ycEng)
        ycInc = zeros(self.n_yc, size(ycEng, 2));
        for i=1:size(ycEng, 2)
            ycInc(:,i) = self.Nyc.*(ycEng(:,i) - self.yc0);
        end
    end
    % Transforms a scaled and incremental ym into one in engineering units
    function ymEng = ym2Eng(self, ymInc)
        ymEng = zeros(self.n_ym, size(ymInc, 2));
        for i=1:size(ymInc, 2)
            ymEng(:,i) = (1./self.Nym).*ymInc(:,i) + self.ym0;
        end
    end
    % Transforms an yr in engineering units into one in scaled and incremental units
    function ymInc = ym2Inc(self, ymEng)
        ymInc = zeros(self.n_ym, size(ymEng, 2));
        for i=1:size(ymEng, 2)
            ymInc(:,i) = self.Nym.*(ymEng(:,i) - self.ym0);
        end
    end
    
    function x_p = x_plus(self, x, u, dm, du)
        % Compute x+ = Ax + Buu*u + Bdu*du + Bdm*dm for the given values of x, u, dm and du
        if nargin == 3
            dm = zeros(self.n_dm, 1);
            du = zeros(self.n_du, 1);
        elseif nargin == 4
            du = zeros(self.n_du, 1);
        end
        % Calculate x+ = Ax + Buu*u + Bdu*du + Bdm*dm
        x_p = self.A*x + self.Bu*u;
        if ~isempty(dm) && ~isempty(self.dm)
            x_p = x_p + self.Bdm*dm;
        end
        if ~isempty(du) && ~isempty(self.du)
            x_p = x_p + self.Bdu*du;
        end
    end
    
    end
    
end

