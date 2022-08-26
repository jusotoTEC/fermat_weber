function Rj=R_j(B,w,j)
  % Operador Rj in Section 7.1 in paper "Weiszfeld’s method: Old and new
  % results"
  % Sintaxis: Rj=R_j(B,w,j)
  % Input:  B = [b_1 b_2 ... b_p] is a m x p matrix
  %         w = [w_1 w_2 ... w_p] is a vector of size p
  %         j = integer number between 1 and p
  % Output: Rj = vector of size m
  % Reference: - Amir Beck and Shoham Sabach. Weiszfeld’s method: Old and new
  %            results. Journal of Optimization Theory and Applications,
  %            164(1):1–40, 2015.
  p=size(B,2);
  Rj=0;
  for i=1:p
    if i~=j
      Rj=Rj+(w(i)/norm(B(:,j)-B(:,i)))*(B(:,j)-B(:,i));
    end
  end
end

