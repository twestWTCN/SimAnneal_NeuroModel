function precomputeBetaSig(R,m)
parfor i = 1:25
    [t,tx_beta,px_beta] = makeBetaSignal(R.model.t_in);
    m_in = m;
    m_in.t = t;
    m_in.tx_beta = tx_beta;
    m_in.px_beta = px_beta;
    saveMkPath([R.rootn '\Inputs\betaSig\betaSignal_' num2str(i)],m_in)
end