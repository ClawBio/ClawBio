import requests
from typing import List, Dict, Any, Optional

# Configuration
BASE_URL = "https://hapi.fhir.org/baseR4"
PATIENT_BATCH_SIZE = 100  # Number of patient IDs per bulk query

def _calculate_age(birth_date_str: str) -> Optional[int]:
    """Return age in years from a FHIR birthDate string (YYYY, YYYY-MM, or YYYY-MM-DD)."""
    if not birth_date_str:
        return None
    try:
        parts = birth_date_str.split("-")
        from datetime import date
        year = int(parts[0])
        month = int(parts[1]) if len(parts) > 1 else 1
        day = int(parts[2]) if len(parts) > 2 else 1
        today = date.today()
        age = today.year - year - ((today.month, today.day) < (month, day))
        return age
    except (ValueError, IndexError):
        return None


def get_all_patients_over_30(base_url: str) -> List[Dict[str, Any]]:
    """
    Retrieve patients older than 30 years from the FHIR server.
    Returns a list of patient dicts with 'id', 'name', and 'age'.
    """
    patients = []
    url = f"{base_url}/Patient"
    params = {"_count": 100}

    while url:
        if url == f"{base_url}/Patient":
            response = requests.get(url, params=params)
        else:
            response = requests.get(url)

        if response.status_code != 200:
            print(f"Error fetching patients: {response.status_code}")
            break

        bundle = response.json()
        for entry in bundle.get("entry", []):
            patient = entry["resource"]
            patient_id = patient["id"]

            age = _calculate_age(patient.get("birthDate", ""))
            if age is None or age <= 30:
                continue

            name = "Unknown"
            if "name" in patient and patient["name"]:
                name_parts = patient["name"][0]
                family = name_parts.get("family", "")
                given = " ".join(name_parts.get("given", []))
                name = f"{family}, {given}".strip(", ")

            patients.append({
                "id": patient_id,
                "name": name,
                "age": age,
            })

        next_link = None
        for link in bundle.get("link", []):
            if link.get("relation") == "next":
                next_link = link.get("url")
                break
        url = next_link

    return patients

def get_resources_for_patients(base_url: str, resource_type: str, patient_ids: List[str], batch_size: int = 100) -> List[Dict]:
    """
    Fetch all resources of a given type (e.g., 'Condition') for a list of patient IDs.
    Handles pagination and returns a flat list of entry resources.
    """
    all_entries = []

    # Split patient IDs into batches to avoid overly long parameters
    for i in range(0, len(patient_ids), batch_size):
        batch_ids = patient_ids[i:i+batch_size]
        patient_refs = [f"Patient/{pid}" for pid in batch_ids]
        patient_param = ",".join(patient_refs)

        # Use POST search to avoid URL length limits
        url = f"{base_url}/{resource_type}/_search"
        data = {"patient": patient_param, "_count": 100}
        headers = {"Content-Type": "application/x-www-form-urlencoded"}

        next_url = url
        while next_url:
            if next_url == url:
                # First request: POST with data
                response = requests.post(next_url, data=data, headers=headers)
            else:
                # Subsequent pages: GET from the next link
                response = requests.get(next_url)

            if response.status_code != 200:
                print(f"Error fetching {resource_type}: {response.status_code}")
                break

            bundle = response.json()
            all_entries.extend(bundle.get("entry", []))

            # Find next page link
            next_link = None
            for link in bundle.get("link", []):
                if link.get("relation") == "next":
                    next_link = link.get("url")
                    break
            next_url = next_link

    return all_entries

def count_resources_per_patient(entries: List[Dict]) -> Dict[str, int]:
    """
    Given a list of bundle entries (each containing a resource), count how many
    resources reference each patient via the 'subject.reference' field.
    """
    counts = {}
    for entry in entries:
        resource = entry.get("resource", {})
        subject = resource.get("subject", {})
        ref = subject.get("reference", "")
        if not ref:
            continue

        # Extract patient ID from reference (assumes format "Patient/123" or full URL)
        if ref.startswith("Patient/"):
            patient_id = ref[8:]
        elif "/Patient/" in ref:
            patient_id = ref.split("/Patient/")[-1]
        else:
            continue  # cannot parse

        counts[patient_id] = counts.get(patient_id, 0) + 1
    return counts

def filter_patients(patients: List[Dict], cond_counts: Dict, req_counts: Dict, stat_counts: Dict,
                    cond_thresh: int = 2, req_thresh: int = 2, stat_thresh: int = 2) -> List[Dict]:
    """
    Return patients that exceed all three thresholds, adding the counts to each patient dict.
    """
    result = []
    for p in patients:
        pid = p["id"]
        c = cond_counts.get(pid, 0)
        r = req_counts.get(pid, 0)
        s = stat_counts.get(pid, 0)
        if c > cond_thresh and r > req_thresh and s > stat_thresh:
            p["condition_count"] = c
            p["medication_request_count"] = r
            p["medication_statement_count"] = s
            result.append(p)
    return result

def main():
    print("Step 1: Retrieving all patients older than 30...")
    patients = get_all_patients_over_30(BASE_URL)
    print(f"Found {len(patients)} patients older than 30.\n")

    if not patients:
        print("No patients to process.")
        return

    patient_ids = [p["id"] for p in patients]

    # Step 2: Fetch clinical data
    print("Step 2: Fetching Condition resources...")
    conditions = get_resources_for_patients(BASE_URL, "Condition", patient_ids, batch_size=PATIENT_BATCH_SIZE)
    print(f"Retrieved {len(conditions)} Condition entries.\n")

    print("Step 2: Fetching MedicationRequest resources...")
    med_requests = get_resources_for_patients(BASE_URL, "MedicationRequest", patient_ids, batch_size=PATIENT_BATCH_SIZE)
    print(f"Retrieved {len(med_requests)} MedicationRequest entries.\n")

    print("Step 2: Fetching MedicationStatement resources...")
    med_statements = get_resources_for_patients(BASE_URL, "MedicationStatement", patient_ids, batch_size=PATIENT_BATCH_SIZE)
    print(f"Retrieved {len(med_statements)} MedicationStatement entries.\n")

    # Step 3: Count resources per patient
    cond_counts = count_resources_per_patient(conditions)
    req_counts = count_resources_per_patient(med_requests)
    stat_counts = count_resources_per_patient(med_statements)

    # Step 4: Apply filter
    filtered = filter_patients(patients, cond_counts, req_counts, stat_counts)

    print(f"\nPatients meeting criteria: {len(filtered)}")
    for p in filtered:
        print(f"ID: {p['id']} | Name: {p['name']} | Age: {p.get('age', '?')} | "
              f"Conditions: {p['condition_count']}, "
              f"MedRequests: {p['medication_request_count']}, "
              f"MedStatements: {p['medication_statement_count']}")

if __name__ == "__main__":
    main()