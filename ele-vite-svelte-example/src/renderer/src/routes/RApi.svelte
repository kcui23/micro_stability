<script>
    let message = '';
    let result = '';
  
    const fetchEcho = async () => {
      const response = await fetch(`http://localhost:8000/echo?msg=${message}`);
      const data = await response.json();
      result = data.msg;
    };
  
    const fetchSum = async () => {
      const response = await fetch(`http://localhost:8000/sum?a=3&b=4`);
      const data = await response.json();
      result = `The sum is: ${data}`;
    };
  </script>
  
  <div>
    <h1>R API Integration</h1>
    <input bind:value={message} placeholder="Enter a message" />
    <button on:click={fetchEcho}>Send to R</button>
    <button on:click={fetchSum}>Fetch Sum</button>
    <p>{result}</p>
  </div>
  