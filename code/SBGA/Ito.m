
function [t, STATES, ALGEBRAIC] = Ito(X, holding_p, holding_t, P1, P1_t, Ek)
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [t, STATES, ALGEBRAIC] = solveModel(X, holding_p, holding_t, P1, P1_t, Ek);
end

function [t, STATES, ALGEBRAIC] = solveModel(X, holding_p, holding_t, P1, P1_t, Ek)
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    
    % Initialise constants and state variables
    INIT_STATES = initConsts();

    % Set timespan to solve over 
    tspan = [0, P1_t];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [t, STATES] = ode15s(@(t, STATES)computeRates(X, t, STATES, holding_p, holding_t, P1, P1_t, Ek), tspan, INIT_STATES, options);
    
    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(X, t, STATES, holding_p, holding_t, P1, P1_t, Ek);
    ALGEBRAIC = computeAlgebraic(X, ALGEBRAIC, STATES, t, holding_p, holding_t, P1, P1_t, Ek);
end

function [RATES, ALGEBRAIC] = computeRates(X, t, STATES, holding_p, holding_t, P1, P1_t, Ek)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    
    % externally applied voltage (voltage clamp)
    ALGEBRAIC(:,6) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t), t);
    
    % Ito
    % A71; alpha_a
    ALGEBRAIC(:,1) =  0.180640.*exp( 0.0357700.*(ALGEBRAIC(:,6)+X(1)));
    % A72; beta_a
    ALGEBRAIC(:,2) =  0.395600.*exp(  - 0.0623700.*(ALGEBRAIC(:,6)+X(2)));
    % A73; alpha_i
    ALGEBRAIC(:,3) = ( 0.000152000.*exp( - (ALGEBRAIC(:,6)+X(3))./7.0))./( 0.00670830.*exp( - (ALGEBRAIC(:,6)+33.5)./7.0)+1.00000);
    % A74; beta_i
    ALGEBRAIC(:,4) = ( 0.000950000.*exp((ALGEBRAIC(:,6)+X(4))./X(5)))./( 0.0513350.*exp((ALGEBRAIC(:,6)+X(4))./X(5))+1.00000);
    % A69; ato_f
    RATES(:,1) =  ALGEBRAIC(:,1).*(1.00000 - STATES(:,1)) -  ALGEBRAIC(:,2).*STATES(:,1);
    % A70; ito_f
    RATES(:,2) =  ALGEBRAIC(:,3).*(1.00000 - STATES(:,2)) -  ALGEBRAIC(:,4).*STATES(:,2);
    % A67; I_Kto,f
    ALGEBRAIC(:,5) =  X(6).*power(STATES(:,1), 3.00000).*STATES(:,2).*(ALGEBRAIC(:,6) - Ek);
    
    RATES = RATES';
end

function ALGEBRAIC = computeAlgebraic(X, ALGEBRAIC, STATES, t, holding_p, holding_t, P1, P1_t, Ek)
    ALGEBRAIC(:,6) = arrayfun(@(t) volt_clamp(t, holding_p, holding_t, P1, P1_t), t);

    % A71; alpha_a
    ALGEBRAIC(:,1) =  0.180640.*exp( 0.0357700.*(ALGEBRAIC(:,6)+X(1)));
    % A72; beta_a
    ALGEBRAIC(:,2) =  0.395600.*exp(  - 0.0623700.*(ALGEBRAIC(:,6)+X(2)));
    % A73; alpha_i
    ALGEBRAIC(:,3) = ( 0.000152000.*exp( - (ALGEBRAIC(:,6)+X(3))./7.0))./( 0.00670830.*exp( - (ALGEBRAIC(:,6)+33.5)./7.0)+1.00000);
    % A74; beta_i
    ALGEBRAIC(:,4) = ( 0.000950000.*exp((ALGEBRAIC(:,6)+X(4))./X(5)))./( 0.0513350.*exp((ALGEBRAIC(:,6)+X(4))./X(5))+1.00000);
    % A67; I_Kto,f
    ALGEBRAIC(:,5) =  X(6).*power(STATES(:,1), 3.00000).*STATES(:,2).*(ALGEBRAIC(:,6) - Ek);
end

function VC = volt_clamp(t, holding_p, holding_t, P1, P1_t)
    if t < holding_t
        VC = holding_p;
    elseif (t >= holding_t) && (t <= P1_t) 
        VC = P1;
    end
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % There are a total of 41 entries in each of the rate and state variable arrays.
    % There are a total of 73 entries in the constant variable array.
    algebraicVariableCount = 6;
end

function [STATES] = initConsts()
    STATES = [];

    STATES(:,1) = 0.265563e-2;  % ato_f; Gating variable for transient outward K+ current
    STATES(:,2) = 0.999977;  % ito_f; Gating variable for transient outward K+ current
    % CONSTANTS(:,1) = 0.4067;  % GKtof; Maximum transient outward K+ current conductance(apex):mS/uF

    if (isempty(STATES)), warning('Initial values for states not set'); end
end
