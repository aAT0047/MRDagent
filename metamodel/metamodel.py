import pandas as pd
import numpy as np
import argparse
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import StratifiedKFold

# ------------------------ Argument Parser ------------------------
parser = argparse.ArgumentParser(description="Train CNN with 10-fold CV on feature and label Excel files.")
parser.add_argument('--features', type=str, required=True, help='Path to Excel file with features (X)')
parser.add_argument('--labels', type=str, required=True, help='Path to Excel file with labels (y)')
args = parser.parse_args()

# ------------------------ Data Load ------------------------
df_X = pd.read_excel(args.features)
df_y = pd.read_excel(args.labels)

X_raw = df_X.values.astype(np.float32)
y_raw = df_y.iloc[:, 0].values

unique_classes = np.unique(y_raw)
if len(unique_classes) > 2:
    is_multiclass = True
    num_classes   = len(unique_classes)
    criterion     = nn.CrossEntropyLoss()
    y_encoded     = y_raw.astype(int)
else:
    is_multiclass = False
    num_classes   = 1
    criterion     = nn.BCEWithLogitsLoss()
    y_encoded     = y_raw.astype(np.float32)

class CustomDataset(Dataset):
    def __init__(self, X, y, multiclass):
        self.X = torch.tensor(X, dtype=torch.float32)
        if multiclass:
            self.y = torch.tensor(y, dtype=torch.long)
        else:
            self.y = torch.tensor(y, dtype=torch.float32)
        self.multiclass = multiclass
    def __len__(self):
        return len(self.y)
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

class CNN1D(nn.Module):
    def __init__(self, input_length, num_classes):
        super().__init__()
        self.conv1 = nn.Conv1d(1, 32, kernel_size=3)
        self.pool1 = nn.MaxPool1d(2)
        self.drop1 = nn.Dropout(0.25)
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3)
        self.pool2 = nn.MaxPool1d(2)
        self.drop2 = nn.Dropout(0.25)
        l1 = (input_length - 2)
        l2 = l1 // 2
        l3 = (l2 - 2)
        l4 = l3 // 2
        self.flatten_dim = 64 * l4
        self.fc1 = nn.Linear(self.flatten_dim, 128)
        self.drop3 = nn.Dropout(0.5)
        self.fc2 = nn.Linear(128, num_classes)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = F.relu(self.conv1(x))
        x = self.pool1(x)
        x = self.drop1(x)
        x = F.relu(self.conv2(x))
        x = self.pool2(x)
        x = self.drop2(x)
        x = x.view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        x = self.drop3(x)
        x = self.fc2(x)
        return x

skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
fold_accuracies = []

for fold, (train_idx, test_idx) in enumerate(skf.split(X_raw, y_raw), 1):
    X_train, X_test = X_raw[train_idx], X_raw[test_idx]
    y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]
    train_dataset = CustomDataset(X_train, y_train, is_multiclass)
    test_dataset  = CustomDataset(X_test, y_test, is_multiclass)
    train_loader  = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader   = DataLoader(test_dataset, batch_size=32)
    model = CNN1D(input_length=X_raw.shape[1], num_classes=num_classes).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    model.train()
    for epoch in range(20):
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            outputs = model(X_batch)
            if is_multiclass:
                loss = criterion(outputs, y_batch)
            else:
                outputs = outputs.squeeze(1)
                loss = criterion(outputs, y_batch)
            loss.backward()
            optimizer.step()
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for X_batch, y_batch in test_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            outputs = model(X_batch)
            if is_multiclass:
                preds = outputs.argmax(dim=1)
            else:
                probs = torch.sigmoid(outputs.squeeze(1))
                preds = (probs >= 0.5).long()
            total += y_batch.size(0)
            correct += (preds == y_batch).sum().item()
    acc = correct / total
    print(f"Fold {fold} Accuracy: {acc:.4f}")
    fold_accuracies.append(acc)

print(f"Average 10-fold Accuracy: {np.mean(fold_accuracies):.4f}")
