function [S_cua] = dimHStabilizer(CMy,rho,v,S,mac,CL_ala,xcg,xac_ala,CL_cua,xac_cua)

syms S_cua;

M_le_ala = CMy *  0.5 * rho * v^2 * S * mac;
L_ala = CL_ala * 0.5 * rho * v^2 * S;
d_XacXcg_ala = xcg - xac_ala;

M_le_cua = 0;
L_cua = CL_cua * 0.5 * v^2 * S_cua;
d_XacXcg_cua = xcg - xac_cua;

Moment_ala = M_le_ala + L_ala * d_XacXcg_ala;
Moment_cua = M_le_cua + L_cua * d_XacXcg_cua;

eq = Moment_ala == Moment_cua;

S_cua = solve(eq,S_cua);

end





