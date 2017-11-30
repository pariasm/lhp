load rnlm.vari.table
rnlm = rnlm_novari;

% each row of the table has
% sigma psz wsz whx wht whvt lambda psnr
%
% out of these, psz, wsz and whvt are constant

for s = [10,20,40],
	
	rnlm_s = rnlm(find(rnlm(:,1) == s),:);
	t = rnlm_s(:,[4,5,6,7,8]);
	
	% find best parameters
	[psnr_max,imax] = max(rnlm_s(:,8));
	best_prms = rnlm_s(find(rnlm_s(:,8) > psnr_max - 0.004*s),:);
	best_prms_ave = mean(best_prms);
	best_prms_std =  std(best_prms);
	
	disp('best parameters (sigma, psz, wsz, hx, ht, hv, lambda, psnr)')
	disp(best_prms)
	
	disp('average and std. dev')
	disp(best_prms_ave)
	disp(best_prms_std)
end

if 0

%hxg = [2*s*[0:.02:1],(4*s*s - 2*s)*[.02:.02:1] + 2*s];
%htg = [2*s*[0:.02:1],(4*s*s - 2*s)*[.02:.02:1] + 2*s];
hxg = 2*s*[0:.05:1];
htg = 2*s*[0:.05:1];
hvg = 2*s*[0:.02:1.3];
lg  = [0:.05:1];
[hx,ht,hv,l] = ndgrid(hxg, htg, hvg, lg);

if ~exist('P'),
	disp('computing matrix P ...')
	P = reshape(griddatan(t(:,1:4), t(:,5), [hx(:), ht(:), hv(:), l(:)], 'linear' ), size(hx));
%	P = reshape(griddatan(t(:,1:3), t(:,4), [hx(:), ht(:), l(:)], 'linear' ), size(hx));
else
	disp('P already computed ... omitting computation')
end

% plot PSNR as a function of ht and hx while varying hv. l is set to 0.85,
% which is the l value of the best result.
mm = min(P(:));
MM = max(P(:));
%mm = MM-2;
for i = 2:length(hvg)-1
	imagesc(hxg, htg, P(:,:,i,18),[mm MM]),
	axis xy
	ylabel('h_t')
	xlabel('h_x')
	axis equal
	axis tight
	colorbar
	title(sprintf('hv = %f (l = l_max)',hvg(i)))

	pause,
end

% plot PSNR as a function of ht and hx while varying hl.
mm = min(P(:));
MM = max(P(:));
%mm = MM-2;
for il = 2:length(lg)-1
	subplot(1,2,1)
	imagesc(hxg, htg, max(P(:,:,:,il),[],3),[mm MM]),
	axis xy
	ylabel('h_t')
	xlabel('h_x')
	axis equal
	axis tight
	colorbar
	title(sprintf('lambda = %f (max wrt hv)',lg(il)))

	subplot(1,2,2)
	imagesc(hxg, htg, min(P(:,:,:,il),[],3),[mm MM]),
	axis xy
	ylabel('h_t')
	xlabel('h_x')
	axis equal
	axis tight
	colorbar
	title(sprintf('lambda = %f (max wrt hv)',lg(il)))
	pause,
end
end

%% % best results for sigma = 10
%% b_hx = hxg(18);
%% b_ht = htg(57);
%% b_l  =  lg(18);
%% 
%% % best results for sigma = 20
%% b_ht = 976;
%% b_hx = 304;
%% b_l  = .95;
%% 
%% % best results for sigma = 40
%% b_ht = 2496;
%% b_hx = 1280;
%% b_l  = .95;
%% 
%% disp('best values:');
%% disp([b_hx b_ht b_l]);
%% 
%% % show table values with similar parameters
%% t(:,1:2) = sqrt(t(:,1:2));
%% disp(t(find((abs(t(:,1) - sqrt(b_hx)) < 5) & ...
%%             (abs(t(:,2) - sqrt(b_ht)) < 5) & ...
%%             (abs(t(:,3) - b_l ) < .3)),:))



% best parameters
%                hx   ht    l   psnr
% sigma = 10 |   68  224  .85   32.5
% sigma = 20 |  304  976  1.0   28.3 
% sigma = 40 | 1608 2995  1.0   24.5
