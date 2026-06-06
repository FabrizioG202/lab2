import ankh
import numpy as np
import torch
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from tqdm import tqdm

import src.data_collection
import src.graphics

# collect data
all_positive, all_negative, positive, negative = src.data_collection.collect_data()

model, tokenizer = ankh.load_base_model()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device

# move to cuda
model = model.to(device)  # ty:ignore[invalid-argument-type]


def get_embeddings(sequences: list[str]) -> list[np.ndarray]:
    tokenized_sequences = [list(seq) for seq in sequences]

    outputs = tokenizer(
        tokenized_sequences,
        add_special_tokens=True,
        padding=True,
        is_split_into_words=True,
        return_tensors="pt",
    )  # ty:ignore[call-non-callable]

    outputs = {k: v.to(device) for k, v in outputs.items()}

    with torch.no_grad():
        model_output = model(
            input_ids=outputs["input_ids"],
            attention_mask=outputs["attention_mask"],
        )

        token_embeddings = model_output.last_hidden_state
        attention_mask = outputs["attention_mask"]

        mask = attention_mask.unsqueeze(-1).float()

        sequence_embeddings = (token_embeddings * mask).sum(dim=1) / mask.sum(
            dim=1
        ).clamp(min=1e-9)

    return [emb for emb in sequence_embeddings.cpu().numpy()]


positive_sequences = positive["sequence"].to_list()
positive_embeddings = []


# loop over the positive sequences in batches of 16 and get the embeddings for each batch
for i in tqdm(range(0, len(positive_sequences), 16)):
    batch_sequences = [s[:90] for s in positive_sequences[i : i + 16]]

    batch_embeddings = get_embeddings(batch_sequences)
    positive_embeddings.extend(batch_embeddings)


# now get them for the negative sequences as well
negative_sequences = negative["sequence"].to_list()
negative_embeddings = []

for i in tqdm(range(0, len(negative_sequences), 16)):
    batch_sequences = [s[:90] for s in negative_sequences[i : i + 16]]

    batch_embeddings = get_embeddings(batch_sequences)
    negative_embeddings.extend(batch_embeddings)

# now, train an mlp on the embeddings to predict the label (positive or negative)


# convert the list of embeddings to a numpy array
positive_embeddings = np.array(positive_embeddings)
negative_embeddings = np.array(negative_embeddings)

# create the labels for the positive and negative samples
positive_labels = np.ones(len(positive_embeddings))
negative_labels = np.zeros(len(negative_embeddings))

# concatenate the positive and negative embeddings and labels
X = np.concatenate((positive_embeddings, negative_embeddings), axis=0)
y = np.concatenate((positive_labels, negative_labels), axis=0)

# split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# train an mlp on the training data
mlp = MLPClassifier(hidden_layer_sizes=(100,), max_iter=1000, random_state=42)
mlp = mlp.fit(X_train, y_train)

# create a cnfusion matrix
# from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay


cm = src.graphics.ConfusionMatrix.from_values(
    y_test,
    mlp.predict(X_test),
)

# Describe the confusion matrix metrics
print(cm.describe())

# save the confusion matrix plot
cm.plot().write_image("./report/.imgs/mlp_confusion_matrix.svg")
