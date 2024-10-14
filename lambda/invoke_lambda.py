import boto3
import json

# my lambda is running on 'eu-north-1' region, which is different from my EC2 region
# All params like aws key, aws secret are taken from my AWS CLI Configuration, then region is set to eu-north-1
lambda_client = boto3.client('lambda', region_name='eu-north-1')

# Define the event you want to send
event_payload = {
    "names": ["Alina", "Raul", "Chris"]
}


# Invoke the Lambda function
response = lambda_client.invoke(
    FunctionName='AlinaZFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event_payload),
)

# Read the response
response_payload = json.loads(response['Payload'].read())
expected_status_code = 200
expected_body = [
    "Hello, Alina!",
    "Hello, Raul!",
    "Hello, Chris!"
]

assert response_payload.get('statusCode') == expected_status_code, f"Expected {expected_status_code}, got {response_payload.get('statusCode')}"
assert response_payload.get('body') == expected_body, f"Expected {expected_body}, got {response_payload.get('body')}"

print("Invoke lambda worked out successfully, unit test passed!")
