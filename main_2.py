import ankh
import torch
from tqdm import tqdm

import src.data_collection

# collect data
all_positive, all_negative, positive, negative = src.data_collection.collect_data()

model, tokenizer = ankh.load_base_model()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device

model = model.to(device)  # ty:ignore[invalid-argument-type]
model.eval()


def get_embeddings(sequences: list[str]):
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
        embeddings = model(
            input_ids=outputs["input_ids"],
            attention_mask=outputs["attention_mask"],
        )

    return embeddings


positive_sequences = positive["sequence"].to_list()
positive_embeddings = []


# loop over the positive sequences in batches of 16 and get the embeddings for each batch
for i in tqdm(range(0, len(positive_sequences), 16)):
    batch_sequences = [s[:90] for s in positive_sequences[i : i + 4]]

    batch_embeddings = get_embeddings(batch_sequences)
    # positive_embeddings.append(batch_embeddings)
