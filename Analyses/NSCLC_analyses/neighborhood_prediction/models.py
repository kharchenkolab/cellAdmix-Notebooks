import torch

def get_model(model_name, params):
    if model_name == 'multilabel_mlp':
        model = MultilabelMLP(params['in_dim'], params['out_dim'], depth=params['enc_depth'], width=params['enc_width'], p_drop=params['p_drop'])
    else:
        raise NotImplementedError
    return model

class MultilabelMLP(torch.nn.Module):
    
    def __init__(self, in_dim, out_dim, depth=4, width=128, p_drop=0.0):
        super(MultilabelMLP, self).__init__()
        self.in_dim = in_dim
        self.out_dim = out_dim
        self.depth = depth
        self.width = width
        if depth > 0:
            layers = []
            layers.append(torch.nn.Dropout(p=p_drop))
            layers.append(torch.nn.Linear(in_dim, width))
            layers.append(torch.nn.ReLU(inplace=True))
            for i in range(depth):
                layers.append(torch.nn.Linear(width, width))
                layers.append(torch.nn.ReLU(inplace=True))
            self.feat_extractor = torch.nn.Sequential(*layers)
            self.head = torch.nn.Linear(width, out_dim)
        else:
            layers = [
                torch.nn.Dropout(p=p_drop),
                torch.nn.Identity()
            ]
            self.feat_extractor = torch.nn.Sequential(*layers) 
            self.head = torch.nn.Linear(in_dim, out_dim)
    
    def forward(self, x, return_feats=False):
        feats = self.feat_extractor(x)
        if return_feats:
            return feats
        logits = self.head(feats)
        return logits
