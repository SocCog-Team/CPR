function [pred_js_dir, pred_ecc] = pred_player(onnx_path, rdp_dir, coherence, acc_mean, ecc_mean, noise, use_ts)
    % Prep data
    rdp_dir      = rdp_dir(:);
    coherence    = coherence(:);
    js_dev       = 360 - ((1 + acc_mean) * 180);
    sin_dev      = sind(js_dev);
    cos_dev      = cosd(js_dev);
    window       = 128;
    N            = numel(rdp_dir);
    if use_ts
        onnx_path = onnx_path(1:length(onnx_path)-5) + "_ts" + onnx_path(length(onnx_path)-4:end);
    end
    
    if use_ts
        % create rdp_dir timestamps
        ts_dir = zeros(N,1);
        run    = 0;    
        for k = 2:N
            if rdp_dir(k) ~= rdp_dir(k-1)
                run = 0;
            else
                run = run + 1;
            end
            ts_dir(k) = run;
        end
        
        s        = sind(rdp_dir);
        c        = cosd(rdp_dir);
        firstRow = [s(1), c(1), coherence(1), ts_dir(1)];
        padding  = repmat(firstRow, window - 1, 1);
        X_block  = [padding; [s, c, coherence, ts_dir]];
    else
        s        = sind(rdp_dir);
        c        = cosd(rdp_dir);
        firstRow = [s(1), c(1), coherence(1)];
        padding  = repmat(firstRow, window - 1, 1);
        X_block  = [padding; [s, c, coherence]];
    end

    Xcell = cell(1, N);
    for n = 1:N
        Xwin     = X_block(n : n + window - 1, :).';
        Xcell{n} = Xwin;
    end
    
    if use_ts
        in_dim = 4;
    else
        in_dim = 3;
    end
    meta_dim = 3;
    B        = numel(Xcell);
    T        = window;
    Xbtc     = zeros(B, T, in_dim + meta_dim, 'single');
    for b = 1:B
        input                = single(Xcell{b});
        Xbtc(b, :, 1:in_dim) = permute(input, [2 1]);
        Xbtc(b, :, in_dim+1) = sin_dev;
        Xbtc(b, :, in_dim+2) = cos_dev;
        Xbtc(b, :, in_dim+3) = ecc_mean;
    end
    
    % Predict player
    params     = importONNXFunction(onnx_path, "LSTM");
    Y          = LSTM(Xbtc, params);
    pred_js_dir = zeros(1, length(Y));
    pred_ecc = zeros(1, length(Y));
    
    for i = 1:length(Y)
        pred_js_dir(i) = mod(atan2d(Y(2, i), Y(1, i)), 360);
        pred_ecc(i)     = sqrt(Y(1, i)^2 + Y(2, i)^2);
    end
    
    if noise > 0
        noise_js_dir = randn(size(pred_js_dir)) * (noise * 360.0);
        noise_ecc = randn(size(pred_ecc)) * noise;

        pred_js_dir = mod(pred_js_dir + noise_js_dir, 360);
        pred_ecc = mod(pred_ecc + noise_ecc, 1);
    end
end